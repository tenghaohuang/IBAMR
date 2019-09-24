// Filename: JacobianCalculator.h
// Created on June 27, 2019 by David Wells and Jordan Brown
//
// Copyright (c) 2019-2019, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "ibtk/JacobianCalculator.h"
#include "ibtk/namespaces.h"

#include "tbox/Utilities.h"

#include <libmesh/point.h>
#include <libmesh/quadrature.h>

#include <algorithm>

namespace
{
// TODO: does libMesh provide this function somewhere?
int
get_dim(const ElemType elem_type)
{
    switch (elem_type)
    {
    case ElemType::TRI3:
    case ElemType::TRI6:
    case ElemType::QUAD4:
    case ElemType::QUAD8:
    case ElemType::QUAD9:
        return 2;
    case ElemType::TET4:
    case ElemType::TET10:
    case ElemType::HEX8:
    case ElemType::HEX27:
        return 3;
    default:
        TBOX_ERROR("unimplemented element type");
    }
    // bogus return to placate compilers
    return 3;
}
} // namespace

JacobianCalculator::JacobianCalculator(const JacobianCalculator::key_type quad_key) : quad_key(quad_key)
{
    const ElemType elem_type = std::get<0>(quad_key);
    const QuadratureType quad_type = std::get<1>(quad_key);
    const Order order = std::get<2>(quad_key);

    const int dim = get_dim(elem_type);

    std::unique_ptr<QBase> quad_rule = QBase::build(quad_type, dim, order);
    quad_rule->init(elem_type);
    d_quad_points = quad_rule->get_points();
    d_quad_weights = quad_rule->get_weights();
    d_JxW.resize(d_quad_weights.size());
}

const std::vector<double>&
JacobianCalculator::get_JxW(const Elem*)
{
    TBOX_ERROR("This base class function is not implemented.");

    return d_JxW;
}

const std::vector<double>&
Tri3JacobianCalculator::get_JxW(const Elem* elem)
{
#ifndef NDEBUG
    TBOX_ASSERT(elem->type() == std::get<0>(quad_key));
#endif
    std::copy(d_quad_weights.begin(), d_quad_weights.end(), d_JxW.begin());

    // calculate Jacobians here
    const Point p0 = elem->point(0);
    const Point p1 = elem->point(1);
    const Point p2 = elem->point(2);

    const double Jac_00 = p1(0) - p0(0);
    const double Jac_01 = p2(0) - p0(0);
    const double Jac_10 = p1(1) - p0(1);
    const double Jac_11 = p2(1) - p0(1);

    const double J = Jac_00 * Jac_11 - Jac_01 * Jac_10;

    TBOX_ASSERT(J > 0.0);
    for (double& jxw : d_JxW) jxw *= J;

    return d_JxW;
}

const std::vector<double>&
Quad4JacobianCalculator::get_JxW(const Elem* elem)
{
#ifndef NDEBUG
    TBOX_ASSERT(elem->type() == std::get<0>(quad_key));
#endif
    std::copy(d_quad_weights.begin(), d_quad_weights.end(), d_JxW.begin());

    // calculate constants in Jacobians here
    const Point p0 = elem->point(0);
    const Point p1 = elem->point(1);
    const Point p2 = elem->point(2);
    const Point p3 = elem->point(3);

    const double a_1 = 0.25 * (-p0(0) + p1(0) + p2(0) - p3(0));
    const double b_1 = 0.25 * (-p0(0) - p1(0) + p2(0) + p3(0));
    const double c_1 = 0.25 * (p0(0) - p1(0) + p2(0) - p3(0));
    const double a_2 = 0.25 * (-p0(1) + p1(1) + p2(1) - p3(1));
    const double b_2 = 0.25 * (-p0(1) - p1(1) + p2(1) + p3(1));
    const double c_2 = 0.25 * (p0(1) - p1(1) + p2(1) - p3(1));

    for (unsigned int i = 0; i < d_JxW.size(); i++)
    {
        // calculate Jacobians here
        const double x = d_quad_points[i](0);
        const double y = d_quad_points[i](1);

        const double Jac_00 = a_1 + c_1 * y;
        const double Jac_01 = b_1 + c_1 * x;
        const double Jac_10 = a_2 + c_2 * y;
        const double Jac_11 = b_2 + c_2 * x;

        const double J = Jac_00 * Jac_11 - Jac_01 * Jac_10;

        TBOX_ASSERT(J > 0.0);
        d_JxW[i] *= J;
    }

    return d_JxW;
}

const std::vector<double>&
Tet4JacobianCalculator::get_JxW(const Elem* elem)
{
#ifndef NDEBUG
    TBOX_ASSERT(elem->type() == std::get<0>(quad_key));
#endif
    std::copy(d_quad_weights.begin(), d_quad_weights.end(), d_JxW.begin());
    
    // calculate Jacobians here
    const Point p0 = elem->point(0);
    const Point p1 = elem->point(1);
    const Point p2 = elem->point(2);
    const Point p3 = elem->point(3);
    
    const double Jac_00 = p1(0) - p0(0);
    const double Jac_01 = p2(0) - p0(0);
    const double Jac_02 = p3(0) - p0(0);
    const double Jac_10 = p1(1) - p0(1);
    const double Jac_11 = p2(1) - p0(1);
    const double Jac_12 = p3(1) - p0(1);
    const double Jac_20 = p1(2) - p0(2);
    const double Jac_21 = p2(2) - p0(2);
    const double Jac_22 = p3(2) - p0(2);
    
    const double J = Jac_00 * (Jac_11 * Jac_22 - Jac_12 * Jac_21) - Jac_10 * (Jac_01 * Jac_22 - Jac_02 * Jac_21) +
    Jac_20 * (Jac_01 * Jac_12 - Jac_02 * Jac_11);
    
    TBOX_ASSERT(J > 0.0);
    for (double& jxw : d_JxW) jxw *= J;
    
    return d_JxW;
}

const std::vector<double>&
Hex8JacobianCalculator::get_JxW(const Elem* elem)
{
#ifndef NDEBUG
    TBOX_ASSERT(elem->type() == std::get<0>(quad_key));
#endif
    std::copy(d_quad_weights.begin(), d_quad_weights.end(), d_JxW.begin());
    
    // calculate constants in Jacobians here
    const Point p0 = elem->point(0);
    const Point p1 = elem->point(1);
    const Point p2 = elem->point(2);
    const Point p3 = elem->point(3);
    const Point p4 = elem->point(4);
    const Point p5 = elem->point(5);
    const Point p6 = elem->point(6);
    const Point p7 = elem->point(7);
    
    const double a_00 = 0.125 * (-p0(0) + p1(0) + p2(0) - p3(0) - p4(0) + p5(0) + p6(0) - p7(0));
    const double b_00 = 0.125 * (p0(0) - p1(0) + p2(0) - p3(0) + p4(0) - p5(0) + p6(0) - p7(0));
    const double c_00 = 0.125 * (p0(0) - p1(0) - p2(0) + p3(0) - p4(0) + p5(0) + p6(0) - p7(0));
    const double d_00 = 0.125 * (-p0(0) + p1(0) - p2(0) + p3(0) + p4(0) - p5(0) + p6(0) - p7(0));
    const double a_10 = 0.125 * (-p0(1) + p1(1) + p2(1) - p3(1) - p4(1) + p5(1) + p6(1) - p7(1));
    const double b_10 = 0.125 * (p0(1) - p1(1) + p2(1) - p3(1) + p4(1) - p5(1) + p6(1) - p7(1));
    const double c_10 = 0.125 * (p0(1) - p1(1) - p2(1) + p3(1) - p4(1) + p5(1) + p6(1) - p7(1));
    const double d_10 = 0.125 * (-p0(1) + p1(1) - p2(1) + p3(1) + p4(1) - p5(1) + p6(1) - p7(1));
    const double a_20 = 0.125 * (-p0(2) + p1(2) + p2(2) - p3(2) - p4(2) + p5(2) + p6(2) - p7(2));
    const double b_20 = 0.125 * (p0(2) - p1(2) + p2(2) - p3(2) + p4(2) - p5(2) + p6(2) - p7(2));
    const double c_20 = 0.125 * (p0(2) - p1(2) - p2(2) + p3(2) - p4(2) + p5(2) + p6(2) - p7(2));
    const double d_20 = 0.125 * (-p0(2) + p1(2) - p2(2) + p3(2) + p4(2) - p5(2) + p6(2) - p7(2));
    
    for (unsigned int i = 0; i < d_JxW.size(); i++)
    {
        // calculate Jacobians here
        const double x = d_quad_points[i](0);
        const double y = d_quad_points[i](1);
        const double z = d_quad_points[i](2);
        
        const double Jac_00 = a_1 + c_1 * y;
        const double Jac_01 = b_1 + c_1 * x;
        const double Jac_02 = 0;
        const double Jac_10 = a_2 + c_2 * y;
        const double Jac_11 = b_2 + c_2 * x;
        const double Jac_12 = 0;
        const double Jac_20 = 0;
        const double Jac_21 = 0;
        const double Jac_22 = 0;
        
        const double J = Jac_00 * (Jac_11 * Jac_22 - Jac_12 * Jac_21) - Jac_10 * (Jac_01 * Jac_22 - Jac_02 * Jac_21) +
        Jac_20 * (Jac_01 * Jac_12 - Jac_02 * Jac_11);
        
        TBOX_ASSERT(J > 0.0);
        d_JxW[i] *= J;
    }
    
    return d_JxW;
}
