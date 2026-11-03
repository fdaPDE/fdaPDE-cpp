// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>

#ifndef __FDAPDE_RGCCA_BOOTSTRAP_NUMA_H__
#define __FDAPDE_RGCCA_BOOTSTRAP_NUMA_H__

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <vector>
#ifdef __linux__
#    include <linux/mempolicy.h>
#    include <sched.h>
#    include <sys/syscall.h>
#    include <unistd.h>
#endif

namespace fdapde::rgcca::internals {

/// @brief restores the calling thread's affinity after a bootstrap worker finishes
class BootstrapAffinityScope {
   public:
    /// @brief pins this worker when a valid physical cpu has been selected
    explicit BootstrapAffinityScope(int cpu) {
#ifdef __linux__
        if (cpu < 0 || sched_getaffinity(0, sizeof(previous_), &previous_) != 0) return;
        cpu_set_t selected;
        CPU_ZERO(&selected);
        CPU_SET(cpu, &selected);
        restore_ = sched_setaffinity(0, sizeof(selected), &selected) == 0;
#else
        static_cast<void>(cpu);
#endif
    }
    /// @brief prevents duplicated affinity restoration
    BootstrapAffinityScope(const BootstrapAffinityScope&) = delete;
    /// @brief prevents duplicated affinity restoration
    BootstrapAffinityScope& operator=(const BootstrapAffinityScope&) = delete;
    /// @brief restores the original affinity on normal and exceptional exits
    ~BootstrapAffinityScope() {
#ifdef __linux__
        if (restore_ && sched_setaffinity(0, sizeof(previous_), &previous_) != 0) std::terminate();
#endif
    }
   private:
#ifdef __linux__
    cpu_set_t previous_;
    bool restore_ = false;
#endif
};

/// @brief owns immutable input replicas for one bootstrap selection between deflations
class BootstrapNumaCache {
   public:
    /// @brief prepares one shared replica per additional numa node within the allowed physical cpus
    BootstrapNumaCache(const std::vector<const Matrix*>& sources, int workers, bool enabled) {
#ifdef __linux__
        if (!enabled || workers < 2 || sources.empty()) return;
        cpu_set_t allowed;
        if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0) return;
        std::map<int, std::map<std::string, std::vector<int>>> groups;
        std::set<std::pair<int, int>> cores;
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu) {
            if (!CPU_ISSET(cpu, &allowed)) continue;
            const auto base = std::filesystem::path("/sys/devices/system/cpu") / ("cpu" + std::to_string(cpu));
            const int package = read_integer_(base / "topology/physical_package_id");
            const int core = read_integer_(base / "topology/core_id");
            if (package < 0 || core < 0 || !cores.emplace(package, core).second) continue;
            int node = -1;
            std::error_code error;
            for (const auto& entry : std::filesystem::directory_iterator(base, error)) {
                const auto name = entry.path().filename().string();
                if (
                  name.starts_with("node") && name.size() > 4 &&
                  name.find_first_not_of("0123456789", 4) == std::string::npos)
                    node = std::stoi(name.substr(4));
            }
            if (error || node < 0 || node >= static_cast<int>(sizeof(unsigned long) * 8)) return;
            std::string cache;
            std::ifstream(base / "cache/index3/shared_cpu_list") >> cache;
            groups[node][cache].push_back(cpu);
        }
        if (groups.size() < 2) return;
        std::map<int, std::vector<int>> node_cpus;
        for (const auto& [node, caches] : groups) {
            for (std::size_t i = 0;; ++i) {
                bool added = false;
                for (const auto& [key, cpus] : caches) {
                    if (i < cpus.size()) {
                        node_cpus[node].push_back(cpus[i]);
                        added = true;
                    }
                }
                if (!added) break;
            }
        }
        for (std::size_t i = 0; static_cast<int>(cpus_.size()) < workers; ++i) {
            bool added = false;
            for (const auto& [node, cpus] : node_cpus) {
                if (i >= cpus.size() || static_cast<int>(cpus_.size()) == workers) continue;
                cpus_.push_back(cpus[i]);
                nodes_.push_back(node);
                added = true;
            }
            if (!added) {
                cpus_.clear();
                nodes_.clear();
                return;
            }
        }
        int home = -1;
        if (
          syscall(SYS_get_mempolicy, &home, nullptr, 0UL, sources.front()->data(), MPOL_F_ADDR | MPOL_F_NODE) != 0 ||
          !node_cpus.contains(home)) {
            cpus_.clear();
            nodes_.clear();
            return;
        }
        for (int node : nodes_) {
            if (node == home || copies_.contains(node)) continue;
            const BootstrapAffinityScope affinity(node_cpus[node].front());
            auto& copies = copies_[node];
            for (const auto* source : sources) {
                std::shared_ptr<Matrix> copy(new Matrix(*source), ReplicaDeleter {});
                if (!has_default_policy_(*copy)) {
                    copies_.clear();
                    cpus_.clear();
                    nodes_.clear();
                    bytes_ = 0;
                    return;
                }
                std::get_deleter<ReplicaDeleter>(copy)->bound = true;
                if (!bind_interior_(*copy, node)) {
                    copies_.clear();
                    cpus_.clear();
                    nodes_.clear();
                    bytes_ = 0;
                    return;
                }
                bytes_ += static_cast<std::size_t>(copy->size()) * sizeof(double);
                copies.push_back(std::move(copy));
            }
        }
#else
        static_cast<void>(sources);
        static_cast<void>(workers);
        static_cast<void>(enabled);
#endif
    }
    /// @brief returns the selected physical cpu or leaves the worker unpinned
    int cpu(int worker) const { return cpus_.empty() ? -1 : cpus_.at(worker); }
    /// @brief borrows a shared owner for this worker's immutable input or retains the original source
    std::shared_ptr<const Matrix> source(int worker, int block) const {
        if (nodes_.empty()) return {};
        const auto found = copies_.find(nodes_.at(worker));
        return found == copies_.end() ? std::shared_ptr<const Matrix> {} : found->second.at(block);
    }
    /// @brief reports additional matrix storage owned by this bootstrap selection
    std::size_t bytes() const { return bytes_; }
   private:
#ifdef __linux__
    /// @brief reads a topology identifier or reports unavailable topology
    static int read_integer_(const std::filesystem::path& path) {
        int value = -1;
        std::ifstream(path) >> value;
        return value;
    }
    /// @brief computes the complete pages exclusively owned by a matrix allocation
    static std::pair<std::uintptr_t, std::uintptr_t> interior_(const Matrix& matrix) noexcept {
        const auto page = static_cast<std::uintptr_t>(sysconf(_SC_PAGESIZE));
        const auto begin = reinterpret_cast<std::uintptr_t>(matrix.data());
        const auto end = begin + static_cast<std::size_t>(matrix.size()) * sizeof(double);
        return {(begin + page - 1) / page * page, end / page * page};
    }
    /// @brief declines replication rather than overriding pre-existing memory policies
    static bool has_default_policy_(const Matrix& matrix) noexcept {
        const auto [begin, end] = interior_(matrix);
        const auto page = static_cast<std::uintptr_t>(sysconf(_SC_PAGESIZE));
        for (auto address = begin; address < end; address += page) {
            int mode = -1;
            if (syscall(SYS_get_mempolicy, &mode, nullptr, 0UL, address, MPOL_F_ADDR) != 0 || mode != MPOL_DEFAULT)
                return false;
        }
        return true;
    }
    /// @brief applies or removes a policy only on complete private matrix pages
    static bool bind_interior_(const Matrix& matrix, int node) noexcept {
        const auto [begin, end] = interior_(matrix);
        if (end <= begin) return true;
        const unsigned long mask = node < 0 ? 0 : 1UL << node;
        // Linux consumes maxnode - 1 mask bits; pass the complete word plus one
        long result;
        do {
            result = syscall(
              SYS_mbind, begin, end - begin, node < 0 ? MPOL_DEFAULT : MPOL_BIND, node < 0 ? nullptr : &mask,
              node < 0 ? 0UL : sizeof(mask) * 8 + 1, node < 0 ? 0 : MPOL_MF_MOVE | MPOL_MF_STRICT);
        } while (result != 0 && errno == EINTR);
        return result == 0;
    }
    /// @brief removes node binding before the last matrix owner returns its allocation to malloc
    struct ReplicaDeleter {
        bool bound = false;
        /// @brief retains the allocation on reset failure so later allocations cannot inherit its policy
        void operator()(const Matrix* matrix) const noexcept {
            if (bound && !bind_interior_(*matrix, -1)) {
                std::fprintf(
                  stderr, "RGCCA: cannot reset NUMA replica policy (errno=%d); retaining %zu bytes\n", errno,
                  static_cast<std::size_t>(matrix->size()) * sizeof(double));
                return;
            }
            delete matrix;
        }
    };
#endif
    std::vector<int> cpus_, nodes_;
    std::map<int, std::vector<std::shared_ptr<const Matrix>>> copies_;
    std::size_t bytes_ = 0;
};
}   // namespace fdapde::rgcca::internals
#endif
