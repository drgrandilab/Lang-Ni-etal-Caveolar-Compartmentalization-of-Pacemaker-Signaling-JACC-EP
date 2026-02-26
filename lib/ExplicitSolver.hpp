/*
	ExplicitSolver_HPP

	explicit solvers to ODEs
	1. Forward Euler
	2.

*/

#ifndef ExplicitSolver_HPP
#define ExplicitSolver_HPP

#include <cmath>


template <typename temp_type>
inline temp_type Foward_Euler(temp_type dt, temp_type dy, temp_type y) {
    return (y + (dt) * dy);
}

template <typename temp_type>
inline temp_type Rush_Larsen(temp_type dt, temp_type gate, temp_type inf, temp_type tau) {
    return (inf + (gate - inf) * exp(-(dt) / tau));
}


#endif