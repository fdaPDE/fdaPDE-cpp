#ifndef __EVALUATOR_IMP_H__
#define __EVALUATOR_IMP_H__

template <UInt ORDER>
void Evaluator<ORDER,2,2>::eval(Real* X, Real *Y, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes,2,2> current_element;
	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	UInt search = mesh_.getSearch();						 

	for (int i = 0; i<length; ++i) {
		current_point = Point(X[i],Y[i]);

		if (search == 1) { //use Naive search
			current_element = mesh_.findLocationNaive(current_point);
		} else if (search == 2)  { //use Tree search (default)
			current_element = mesh_.findLocationTree(current_point);
		} else if (search == 3) { //use Walking search
			Element<Nodes,2,2> starting_element;
			starting_element = mesh_.getElement(0);
			current_element = mesh_.findLocationWalking(current_point, starting_element);
			if(current_element.getId() == Identifier::NVAL && redundancy == true) {
				//To avoid problems with non convex mesh
				current_element = mesh_.findLocationNaive(current_point);
			}
		}

		if(current_element.getId() == Identifier::NVAL) {
			isinside[i]=false;
		} else {
			isinside[i]=true;
			for (int j=0; j<(Nodes); ++j) {
				coefficients[j] = coef[current_element[j].getId()];
			}
			// std::cout << "i : " << i << " current element id: " << current_element.getId() << std::endl;
			result[i] = evaluate_point<Nodes,2,2>(current_element, current_point, coefficients);
		}
	} //end of for loop
}

template <UInt ORDER>
void Evaluator<ORDER,2,2>::evalWithInfo(Real* X, Real *Y, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside, const std::vector<UInt> & element_id, Real **barycenters)
{

	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes,2,2> current_element;
	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	Eigen::Matrix<Real,Nodes,1> bary_coeff;

	for (int i = 0; i<length; ++i) {
		current_point = Point(X[i],Y[i]);
		current_element = mesh_.getElement(element_id[i]);

		if(current_element.getId() == Identifier::NVAL) {
			isinside[i]=false;
		} else {
			isinside[i]=true;
			for (int j=0; j<Nodes; ++j) {
				coefficients[j] = coef[current_element[j].getId()];
				bary_coeff[j]=barycenters[i][j];
			}

			if (Nodes == 3) {
				result[i] = coefficients.dot(bary_coeff);
			} else if (Nodes == 6) {
				result[i] = coefficients[0]*(2*bary_coeff[0]*bary_coeff[0] - bary_coeff[0]) +
				            coefficients[1]*(2*bary_coeff[1]*bary_coeff[1] - bary_coeff[1]) +
				            coefficients[2]*(2*bary_coeff[2]*bary_coeff[2] - bary_coeff[2]) +
				            coefficients[3]*(4*bary_coeff[1]* bary_coeff[2]) +
				            coefficients[4]*(4*bary_coeff[2]* bary_coeff[0]) +
				            coefficients[5]*(4*bary_coeff[0]* bary_coeff[1]);
			}
		}
	} //end of for loop
}



template <UInt ORDER>
void Evaluator<ORDER,2,3>::eval(Real* X, Real *Y,  Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes,2,3> current_element;
	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	UInt search = mesh_.getSearch();

	for (int i = 0; i<length; ++i) {
		current_point = Point(X[i],Y[i],Z[i]);

		if (search == 1) { //use Naive search
			current_element = mesh_.findLocationNaive(current_point);
		} else if (search == 2)  { //use Tree search (default)
			current_element = mesh_.findLocationTree(current_point);
		}

		if(current_element.getId() == Identifier::NVAL) {
			isinside[i]=false;
		} else {
			isinside[i]=true;
			for (int j=0; j<(Nodes); ++j) {
				coefficients[j] = coef[current_element[j].getId()];
			}
			// std::cout << "i : " << i << " current element id: " << current_element.getId() << std::endl;
			result[i] = evaluate_point<Nodes,2,3>(current_element, current_point, coefficients);
		}
	} //end of for loop

}

template <UInt ORDER>
void Evaluator<ORDER,2,3>::evalWithInfo(Real* X, Real *Y, Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside, const std::vector<UInt> & element_id, Real **barycenters)
{

	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes,2,3> current_element;
	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	Eigen::Matrix<Real,Nodes,1> bary_coeff;

	for (int i = 0; i<length; ++i) {
		current_point = Point(X[i],Y[i],Z[i]);
		current_element = mesh_.getElement(element_id[i]);

		if(current_element.getId() == Identifier::NVAL) {
			isinside[i]=false;
		} else {
			isinside[i]=true;
			for (int j=0; j<Nodes; ++j) {
				coefficients[j] = coef[current_element[j].getId()];
				bary_coeff[j]=barycenters[i][j];
			}

			if (Nodes == 3) {
				result[i] = coefficients.dot(bary_coeff);
			} else if (Nodes == 6) {
				result[i] = coefficients[0]*(2*bary_coeff[0]*bary_coeff[0] - bary_coeff[0]) +
				            coefficients[1]*(2*bary_coeff[1]*bary_coeff[1] - bary_coeff[1]) +
				            coefficients[2]*(2*bary_coeff[2]*bary_coeff[2] - bary_coeff[2]) +
				            coefficients[3]*(4*bary_coeff[1]*bary_coeff[2]) +
				            coefficients[4]*(4*bary_coeff[2]*bary_coeff[0]) +
				            coefficients[5]*(4*bary_coeff[0]*bary_coeff[1]);
			}

		}
	} //end of for loop
}

template <UInt ORDER>
void Evaluator<ORDER,3,3>::eval(Real* X, Real *Y,  Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside)
{

	constexpr UInt Nodes = 6*ORDER-2;
	Element<Nodes,3,3> current_element;
	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	UInt search = mesh_.getSearch();


	for (int i = 0; i<length; ++i) {
		current_point = Point(X[i],Y[i],Z[i]);

		if (search == 1) { //use Naive search
			current_element = mesh_.findLocationNaive(current_point);
		} else if (search == 2)  { //use Tree search (default)
			current_element = mesh_.findLocationTree(current_point);
		}


		if(current_element.getId() == Identifier::NVAL) {
			isinside[i]=false;
		} else {
			isinside[i]=true;
			for (int j=0; j<(Nodes); ++j) {
				coefficients[j] = coef[current_element[j].getId()];
			}
			// std::cout << "i : " << i << " current element id: " << current_element.getId() << std::endl;
			result[i] = evaluate_point<Nodes,3,3>(current_element, current_point, coefficients);
		}
	} //end of for loop

}

template <UInt ORDER>
void Evaluator<ORDER,3,3>::evalWithInfo(Real* X, Real *Y, Real *Z, UInt length, const Real *coef, bool redundancy, Real* result, std::vector<bool>& isinside, const std::vector<UInt> & element_id, Real **barycenters)
{

	constexpr UInt Nodes = 6*ORDER-2;
	Element<Nodes,3,3> current_element;
	Point current_point;
	Eigen::Matrix<Real,Nodes,1> coefficients;
	Eigen::Matrix<Real,Nodes,1> bary_coeff;

	for (int i = 0; i<length; ++i) {
		current_point = Point(X[i],Y[i],Z[i]);
		current_element = mesh_.getElement(element_id[i]);

		if(current_element.getId() == Identifier::NVAL) {
			isinside[i]=false;
		} else {
			isinside[i]=true;
			for (int j=0; j<Nodes; ++j) {
				coefficients[j] = coef[current_element[j].getId()];
				bary_coeff[j]=barycenters[i][j];
			}

			result[i] = coefficients.dot(bary_coeff);
		}
	} //end of for loop
}


template <UInt ORDER>
void Evaluator<ORDER, 2, 2>::integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result)
{
	std::vector<Real> Delta(nRegions);
	std::vector<Real> integral(nRegions);
	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes, 2, 2> current_element;

	for (int region=0; region<nRegions; region++)
	{
		Delta[region]=0;
		integral[region]=0;
		for (int elem=0; elem<nElements; elem++)
		{
			if (incidenceMatrix[region][elem]==1) //elem is in region
			{
				current_element = mesh_.getElement(elem);
				Real measure = mesh_.elementMeasure(elem);
				Delta[region] += measure;
				Real s = 0;
				for (int node = ORDER==1 ? 0 : 3; node<Nodes; node++)
				{
					s+=coef[current_element[node].getId()];
				}
				integral[region] += measure*s/(2+1);
			}
		}
		result[region]=integral[region]/Delta[region];
	}
}


template <UInt ORDER>
void Evaluator<ORDER, 2, 3>::integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result)
{
	std::vector<Real> Delta(nRegions);
	std::vector<Real> integral(nRegions);
	constexpr UInt Nodes = 3*ORDER;
	Element<Nodes, 2, 3> current_element;

	for (int region=0; region<nRegions; region++)
	{
		Delta[region]=0;
		integral[region]=0;
		for (int elem=0; elem<nElements; elem++)
		{
			if (incidenceMatrix[region][elem]==1) //elem is in region
			{

				current_element = mesh_.getElement(elem);
				Real measure = mesh_.elementMeasure(elem);
				Delta[region] += measure;

				Real s = 0;
				for (int node = ORDER==1 ? 0 : 3; node<Nodes; node++)
				{
					s+=coef[current_element[node].getId()];
				}
				integral[region] += measure*s/(2+1);

			}
		}
		result[region]=integral[region]/Delta[region];

	}
}


template <UInt ORDER>
void Evaluator<ORDER, 3, 3>::integrate(UInt** incidenceMatrix, UInt nRegions, UInt nElements, const Real *coef, Real* result)
{
	std::vector<Real> Delta(nRegions);
	std::vector<Real> integral(nRegions);
	constexpr UInt Nodes = 6*ORDER-2;
	Element<Nodes, 3, 3> current_element;

	for (int region=0; region<nRegions; region++)
	{
		Delta[region]=0;
		integral[region]=0;
		for (int elem=0; elem<nElements; elem++)
		{
			if (incidenceMatrix[region][elem]==1) //elem is in region
			{
				current_element = mesh_.getElement(elem);
				Real measure = mesh_.elementMeasure(elem);
				Delta[region] += measure;

				// THIS IS ONLY FOR ORDER==1
				Real s = 0;
				for (int node=0; node<Nodes; node++)
				{
					s+=coef[current_element[node].getId()];
				}
				integral[region] += measure*s/(3+1);

			}
		}
		result[region]=integral[region]/Delta[region];

	}
}


#endif
