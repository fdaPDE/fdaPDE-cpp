#ifndef __EVALUATOR_IMP_H__
#define __EVALUATOR_IMP_H__

template<UInt ORDER, UInt mydim, UInt ndim>
template<bool isManifold>
typename std::enable_if<!isManifold,void>::type
Evaluator<ORDER, mydim, ndim>::eval(const RNumericMatrix &locations, const RNumericMatrix &coef, bool redundancy,
                                        RNumericMatrix &result, std::vector<bool> &isinside){
    UInt length= locations.nrows();
    constexpr UInt Nodes = how_many_nodes(ORDER,mydim);
    Element<Nodes,mydim,ndim> current_element;
    Point<ndim> current_point;
    Eigen::Matrix<Real,Nodes,1> coefficients;
    UInt search = mesh_.getSearch();

    for (UInt i = 0; i<length; ++i) {
        std::array<Real, ndim> coords;
        for(UInt n=0; n<ndim; ++n)
            coords[n] = locations(i,n);

        current_point = Point<ndim>(coords);
        current_element = mesh_.findLocation(current_point);
        if(search==3 && current_element.getId() == Identifier::NVAL && redundancy == true) {
            //To avoid problems with non convex mesh
            current_element = mesh_.findLocationNaive(current_point);
        }

        if(current_element.getId() == Identifier::NVAL) {
            isinside[i]=false;
        } else {
            isinside[i]=true;
            for (int j=0; j<(Nodes); ++j) {
                coefficients[j] = coef[current_element[j].getId()];
            }
            // std::cout << "i : " << i << " current element id: " << current_element.getId() << std::endl;
            result[i] = current_element.evaluate_point(current_point, coefficients);
        }
    } //end of for loop
}

template<UInt ORDER, UInt mydim, UInt ndim>
template<bool isManifold>
typename std::enable_if<isManifold,void>::type
Evaluator<ORDER, mydim, ndim>::eval(const RNumericMatrix &locations, const RNumericMatrix &coef, bool redundancy,
                                        RNumericMatrix &result, std::vector<bool> &isinside) {
    UInt length = locations.nrows();
    constexpr
    UInt Nodes = how_many_nodes(ORDER, mydim);
    Element<Nodes, mydim, ndim> current_element;
    Point<ndim> current_point;
    Eigen::Matrix<Real, Nodes, 1> coefficients;
    UInt search = mesh_.getSearch();

    for (UInt i = 0; i < length; ++i) {
        std::array <Real, ndim> coords;
        for (UInt n = 0; n < ndim; ++n)
            coords[n] = locations(i, n);

        current_point = Point<ndim>(coords);
        current_element = mesh_.findLocation(current_point);

        if (current_element.getId() == Identifier::NVAL) {
            isinside[i] = false;
        } else {
            isinside[i] = true;
            for (int j = 0; j < (Nodes); ++j) {
                coefficients[j] = coef[current_element[j].getId()];
            }
            // std::cout << "i : " << i << " current element id: " << current_element.getId() << std::endl;
            result[i] = current_element.evaluate_point(current_point, coefficients);
        }
    } //end of for loop
}

template<UInt ORDER,UInt mydim, UInt ndim>
void Evaluator<ORDER, mydim, ndim>::evalWithInfo(const RNumericMatrix &locations, const RNumericMatrix &coef,
                                                     bool redundancy, RNumericMatrix &result,
                                                     std::vector<bool> &isinside, const RIntegerMatrix &element_id,
                                                     const RNumericMatrix &barycenters) {
    UInt length = locations.nrows();
    constexpr UInt Nodes = how_many_nodes(ORDER,mydim);
    Element<Nodes,mydim,ndim> current_element;
    Point<ndim> current_point;
    Eigen::Matrix<Real,Nodes,1> coefficients;
    UInt search = mesh_.getSearch();

    for (int i = 0; i<length; ++i) {
        std::array<Real, ndim> coords;
        for(UInt n=0; n<ndim;++n)
            coords[n] = locations(i,n);

        current_point = Point<ndim>(coords);
        current_element = mesh_.getElement(element_id[i]);

        if(current_element.getId() == Identifier::NVAL) {
            isinside[i]=false;
        } else {
            isinside[i]=true;
            for (int j=0; j<Nodes; ++j) {
                coefficients[j] = coef[current_element[j].getId()];
            }

            result[i] = current_element.evaluate_point(current_point, coefficients);
        }
    } //end of for loop

}

template<UInt ORDER, UInt mydim, UInt ndim>
void Evaluator<ORDER, mydim, ndim>::integrate(const RIntegerMatrix &incidenceMatrix, const RNumericMatrix &coef,
                                                  RNumericMatrix &result) {
    UInt nRegions = incidenceMatrix.nrows();
    UInt nElements = incidenceMatrix.ncols();
    std::vector<Real> Delta(nRegions);
    std::vector<Real> integral(nRegions);
    constexpr UInt Nodes = how_many_nodes(ORDER,mydim);
    Element<Nodes, mydim, ndim> current_element;
    Eigen::Matrix<Real,Nodes,1> coefficients;

    for (UInt region=0; region<nRegions; ++region)
    {
        for (UInt elem=0; elem<nElements; ++elem)
        {
            if (incidenceMatrix(region,elem)==1) //elem is in region
            {
                current_element = mesh_.getElement(elem);
                for (UInt i=0; i<Nodes; ++i)
                    coefficients[i]=coef[current_element[i].getId()];
                Delta[region] += current_element.getMeasure();
                integral[region] += current_element.integrate(coefficients);
            }
        }
        result[region]=integral[region]/Delta[region];
    }
}


#endif
