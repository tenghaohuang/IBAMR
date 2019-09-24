// Copyright (c) 2019, Boyce Griffith
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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>

#include <SAMRAI_config.h>

// Headers for basic libMesh objects
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>

// Headers for application-specific algorithm/data structure objects
#include <ibtk/AppInitializer.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/JacobianCalculator.h>

// Set up application namespace declarations
#include <ibamr/app_namespaces.h>

#include <boost/multi_array.hpp>

// Verify that JacobianCalc and descendants output the same values as libMesh::FEMap.

int
main(int argc, char** argv)
{
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    {
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");

        Pointer<Database> input_db = app_initializer->getInputDatabase();
        const double radius = 1.0;
        const unsigned int n_refinements = 1;

        {
            plog << "Test 1: TRI3" << std::endl;
            ReplicatedMesh mesh(init.comm(), 2);
            // MeshTools::Generation::build_sphere(mesh, radius, n_refinements, TRI3);
            MeshTools::Generation::build_square(mesh, 10, 10, 0.0, 1.0, 0.0, 1.0, TRI3);
            std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order> key(TRI3, QGAUSS, FOURTH);
            // libMesh::MeshRefinement mesh_refinement(mesh);
            // mesh_refinement.uniformly_refine(4);
            Tri3JacobianCalculator jac_calc(key);
            const std::vector<double> &JxW = jac_calc.get_JxW(*mesh.active_local_elements_begin());
            for (const double jxw : JxW)
                plog << std::setprecision(12) << jxw << '\n';

            double volume = 0;
            for (auto elem_iter = mesh.active_local_elements_begin();
                 elem_iter != mesh.active_local_elements_end(); ++elem_iter)
            {
                const std::vector<double> &JxW = jac_calc.get_JxW(*elem_iter);
                volume += std::accumulate(JxW.begin(), JxW.end(), 0.0);
            }
            plog << "volume is " << volume << '\n';
        }

        {
            plog << "Test 2: QUAD4" << std::endl;
            ReplicatedMesh mesh(init.comm(), 2);
            // MeshTools::Generation::build_sphere(mesh, radius, n_refinements, QUAD4);
            MeshTools::Generation::build_square(mesh, 10, 10, 0.0, 1.0, 0.0, 1.0, QUAD4);
            std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order> key(QUAD4, QGAUSS, FOURTH);
            // libMesh::MeshRefinement mesh_refinement(mesh);
            // mesh_refinement.uniformly_refine(4);
            Quad4JacobianCalculator jac_calc(key);
            const std::vector<double> &JxW = jac_calc.get_JxW(*mesh.active_local_elements_begin());
            for (const double jxw : JxW)
                plog << std::setprecision(12) << jxw << '\n';

            double volume = 0;
            for (auto elem_iter = mesh.active_local_elements_begin();
                 elem_iter != mesh.active_local_elements_end(); ++elem_iter)
            {
                const std::vector<double> &JxW = jac_calc.get_JxW(*elem_iter);
                volume += std::accumulate(JxW.begin(), JxW.end(), 0.0);
            }
            plog << "volume is " << volume << '\n';
        }
        
        {
            plog << "Test 3: TET4" << std::endl;
            ReplicatedMesh mesh(init.comm(), 2);
            MeshTools::Generation::build_cube(mesh, 10, 10, 10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, TET4);
            std::tuple<libMesh::ElemType, libMesh::QuadratureType, libMesh::Order> key(TET4, QGAUSS, FOURTH);
            Tet4JacobianCalculator jac_calc(key);
            const std::vector<double> &JxW = jac_calc.get_JxW(*mesh.active_local_elements_begin());
            for (const double jxw : JxW)
                plog << std::setprecision(12) << jxw << '\n';
            
            double volume = 0;
            for (auto elem_iter = mesh.active_local_elements_begin();
                 elem_iter != mesh.active_local_elements_end(); ++elem_iter)
            {
                const std::vector<double> &JxW = jac_calc.get_JxW(*elem_iter);
                volume += std::accumulate(JxW.begin(), JxW.end(), 0.0);
            }
            plog << "volume is " << volume << '\n';
        }
#if 0
        {
            plog << "Test 2: tri6" << std::endl;
            ReplicatedMesh mesh(init.comm(), 2);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, TRI6);
        }

        {
            plog << std::endl << "Test 3: quad4" << std::endl;
            ReplicatedMesh mesh(init.comm(), 2);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, QUAD4);
        }

        {
            plog << std::endl << "Test 4: quad9" << std::endl;
            ReplicatedMesh mesh(init.comm(), 2);
            MeshTools::Generation::build_sphere(mesh, radius, n_refinements, QUAD9);
        }
#endif
    }

    SAMRAIManager::shutdown();
} // main
