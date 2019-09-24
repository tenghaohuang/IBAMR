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

#ifndef included_IBTK_JacobianCalculator
#define included_IBTK_JacobianCalculator

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/enum_order.h>
#include <libmesh/enum_quadrature_type.h>
#include <libmesh/point.h>

#include <tuple>
#include <vector>

class JacobianCalculator
{
public:
    /**
     * Key type. Completely describes (excepting p-refinement) a libMesh
     * quadrature rule.
     */
    using key_type = std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order>;

    /**
     * Constructor.
     */
    JacobianCalculator(const key_type quad_key);

    /**
     * Calculate the JxW values on the given element and return a reference to
     * the result.
     */
    virtual
    const std::vector<double> &
    get_JxW(const libMesh::Elem * elem) = 0;

protected:
    const key_type quad_key;

    std::vector<libMesh::Point> d_quad_points;
    std::vector<double> d_quad_weights;

    std::vector<double> d_JxW;

};

class Tri3JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual
    const std::vector<double> &
    get_JxW(const libMesh::Elem * elem) override;
};

class Quad4JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual
    const std::vector<double> &
    get_JxW(const libMesh::Elem * elem) override;
};

class Tet4JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual
    const std::vector<double> &
    get_JxW(const libMesh::Elem * elem) override;
};

class Hex8JacobianCalculator : public JacobianCalculator
{
public:
    /**
     * Explicitly use the base class' constructor (this class does not require
     * any additional setup).
     */
    using JacobianCalculator::JacobianCalculator;

    virtual
    const std::vector<double> &
    get_JxW(const libMesh::Elem * elem) override;
};

#endif //#ifndef included_IBTK_JacobianCalculator
