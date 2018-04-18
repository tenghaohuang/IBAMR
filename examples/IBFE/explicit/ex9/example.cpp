// Copyright (c) 2002-2014, Boyce Griffith
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
//     /*---- This example is setup for falling sphere using the sharp IBFE ---*/
//
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

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for basic libMesh objects
#include <libmesh/boundary_info.h>
#include <libmesh/boundary_mesh.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_triangle_interface.h>

// Headers for application-specific algorithm/data structure objects
#include <boost/multi_array.hpp>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LEInteractor.h>
#include <ibtk/libmesh_utilities.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

#include <ibtk/LData.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBStandardInitializer.h>


// Set up application namespace declarations
#include <ibamr/app_namespaces.h>


// Elasticity model data.
namespace ModelData
{
static double kappa_s = 1.0e6;
static double eta_s = 0.0;
static double grav_const[3]={0.0, -9.81,0.0};


System* x_new_solid_system, * u_new_solid_system;
System* x_half_solid_system, * u_half_solid_system;
void
tether_force_function(VectorValue<double>& F,
                      const TensorValue<double>& /*FF*/,
                      const libMesh::Point& x_bndry,  // x_bndry gives current   coordinates on the boundary mesh
                      const libMesh::Point& X_bndry,  // X_bndry gives reference coordinates on the boundary mesh
                      Elem* const elem,
                      const vector<const vector<double>*>& var_data,
                      const vector<const vector<VectorValue<double> >*>& /*grad_var_data*/,
                      double /*time*/,
                      void* /*ctx*/)
{
    // tether_force_function() is called on elements of the boundary mesh.  Here
    // we look up the element in the solid mesh that the current boundary
    // element was extracted from.
    const Elem* const interior_parent = elem->interior_parent();

    // We define "arbitrary" velocity and displacement fields on the solid mesh.
    // Here we look up their values.
    std::vector<double> x_solid(NDIM, 0.0), u_solid(NDIM, 0.0);
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        x_solid[d] = x_new_solid_system->point_value(d, X_bndry, interior_parent);
        u_solid[d] = u_new_solid_system->point_value(d, X_bndry, interior_parent);
    }

    // Look up the velocity of the boundary mesh.
    const std::vector<double>& u_bndry = *var_data[0];

    // The tether force is proportional to the mismatch between the positions
    // and velocities.
    for (unsigned int d = 0; d < NDIM; ++d)
    {
        F(d) = kappa_s * (x_solid[d] - x_bndry(d)) + eta_s * (u_solid[d] - u_bndry[d]);
    }
    return;
} // tether_force_function
}
using namespace ModelData;

static ofstream w_new_stream, v_new_stream, x_com_stream;
void postprocess_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                      Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                      Mesh& mesh,
                      EquationSystems* equation_systems,
                      const int iteration_num,
                      const double loop_time,
                      const string& data_dump_dirname);


/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
void calculateGeomQuantitiesOfStructure(double& M_new,  // mass of the body
										TensorValue<double>& I_w_new,  // moment of inertia tensor
										VectorValue<double>& x_com_new,        // new center of the mass
										const double rho_s,
										EquationSystems* solid_equation_systems)                  // mass density of the body (assumed to be uniform)
										//~ libMesh::UniquePtr<EquationSystems> solid_equation_systems)
{
    // Get the structure mesh for codim-0 solid.
    // For now the eqs are setup only for one part but this will be extended
    // to multiple parts
    

 

    // Extract the FE system and DOF map, and setup the FE object.
    //~ System& X_new_system = solid_equation_systems->get_system("position_new");
    //~ x_new_solid_system->solution->localize(*x_new_solid_system->current_local_solution);
    
                  
     System& X_system = solid_equation_systems->get_system("position_new");
     X_system.solution->localize(*X_system.current_local_solution);
     const unsigned int X_sys_num = X_system.number();
      MeshBase& mesh = solid_equation_systems->get_mesh();
                  
                 
                
      //~ MeshBase& mesh = solid_equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
 
    DofMap& X_dof_map = X_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    FEType fe_type = X_dof_map.variable_type(0);
    
    

    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    PetscVector<double>& X_new_petsc = dynamic_cast<PetscVector<double>&>(*X_system.current_local_solution.get());
    X_new_petsc.close();
    Vec X_new_global_vec = X_new_petsc.vec();
    Vec X_new_local_ghost_vec;
    VecGhostGetLocalForm(X_new_global_vec, &X_new_local_ghost_vec);
    double* X_new_local_ghost_soln;
    VecGetArray(X_new_local_ghost_vec, &X_new_local_ghost_soln);
    

   
    
    // 3D identity tensor.
    static const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

    M_new = 0.0;
    I_w_new.zero();
    x_com_new.zero();
    double vol_new = 0.0; //volume of the body
   boost::multi_array<double, 2> X_new_node;



   // double X_qp_new[NDIM], X_qp_current[NDIM], R_qp_current[NDIM], R_qp_new[NDIM];
    VectorValue<double> X_qp_new, R_qp_new;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);

        }
        get_values_for_interpolation(X_new_node, X_new_petsc, X_new_local_ghost_soln, X_dof_indices);


        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp_new, qp, X_new_node, phi);

            x_com_new += X_qp_new * JxW[qp];
            
            vol_new += JxW[qp];

        }
       
    }
    SAMRAI_MPI::sumReduction(&x_com_new(0), NDIM);
    SAMRAI_MPI::sumReduction(&vol_new);

                        
                        
    x_com_new /= vol_new;
    M_new = rho_s * vol_new;
    
    
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);

        }
        get_values_for_interpolation(X_new_node, X_new_petsc, X_new_local_ghost_soln, X_dof_indices);


        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(X_qp_new, qp, X_new_node, phi);

                      
			R_qp_new = X_qp_new - x_com_new;

            // Accumulate the inertia tensor:
            I_w_new += rho_s * ((R_qp_new * R_qp_new) * II - outer_product(R_qp_new, R_qp_new)) * JxW[qp];


        }
    }
    
     SAMRAI_MPI::sumReduction(&I_w_new(0, 0), NDIM * NDIM);

    VecRestoreArray(X_new_local_ghost_vec, &X_new_local_ghost_soln);
    VecGhostRestoreLocalForm(X_new_global_vec, &X_new_local_ghost_vec);
   
    
     X_system.solution->close();


    
    return;
    
} // calculateGeomQuantitiesOfStructure



void calculateFluidForceAndTorque(VectorValue<double>& F,              // net force  acting on the body
								  VectorValue<double>& T,              // net torque acting on the body
								  VectorValue<double> x_com_current,
								  Mesh& mesh,
								  EquationSystems* equation_systems)
{
	
	 const unsigned int dim = mesh.mesh_dimension();
     F.zero();
     T.zero();
     
    System& x_system = equation_systems->get_system(IBFEMethod::COORDS_SYSTEM_NAME);
    System& TAU_system = equation_systems->get_system<System>(IBFEMethod::TAU_SYSTEM_NAME);
    
    NumericVector<double>* TAU_vec = TAU_system.solution.get();
    NumericVector<double>* TAU_ghost_vec = TAU_system.current_local_solution.get();
    TAU_vec->localize(*TAU_ghost_vec);
    DofMap& TAU_dof_map = TAU_system.get_dof_map();
    std::vector<std::vector<unsigned int> > TAU_dof_indices(NDIM);
    AutoPtr<FEBase> fe(FEBase::build(dim, TAU_dof_map.variable_type(0)));
    
    
    NumericVector<double>* x_vec = x_system.solution.get();
    NumericVector<double>* x_ghost_vec = x_system.current_local_solution.get();
    x_vec->localize(*x_ghost_vec);
    const DofMap& dof_map = x_system.get_dof_map();
    std::vector<std::vector<unsigned int> > dof_indices(NDIM);
    
    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const vector<double>& JxW = fe->get_JxW();
    const vector<vector<double> >& phi = fe->get_phi();
    
    
    boost::multi_array<double, 2> x_node, TAU_node;
    VectorValue<double> F_qp, x_qp, W_qp, TAU_qp, R_qp;
        
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;
        fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            dof_map.dof_indices(elem, dof_indices[d], d);
            TAU_dof_map.dof_indices(elem, TAU_dof_indices[d], d);
        }
        get_values_for_interpolation(x_node, *x_ghost_vec, dof_indices);
        get_values_for_interpolation(TAU_node, *TAU_ghost_vec, TAU_dof_indices);
                      
 
        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {
            interpolate(x_qp, qp, x_node, phi);
            interpolate(TAU_qp, qp, TAU_node, phi);


            
            R_qp = x_qp - x_com_current;
            
            F += TAU_qp * JxW[qp];
            T += R_qp.cross(TAU_qp) * JxW[qp];

        }
         
 
    }
    SAMRAI_MPI::sumReduction(&F(0), NDIM);
    SAMRAI_MPI::sumReduction(&T(0), NDIM);


     x_ghost_vec->close();
     TAU_ghost_vec->close();
     
	return;
	

} //calculateFluidForceAndTorque



void calculateGravitationalForce(VectorValue<double>& F_g, //gravitational body force
								 const double rho,        // mass density of the fluid 
								 const double rho_s,      //mass density of the body (assumed to be uniform)
								 EquationSystems* solid_equation_systems)  
{
	
	MeshBase& mesh = solid_equation_systems->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);


    // Extract the FE system and DOF map, and setup the FE object.
    System& X_half_system = solid_equation_systems->get_system("position_half");


    X_half_system.solution->localize(*X_half_system.current_local_solution);
    DofMap& X_half_dof_map = X_half_system.get_dof_map();
    std::vector<std::vector<unsigned int> > X_half_dof_indices(NDIM);
    
    
    FEType fe_type_half = X_half_dof_map.variable_type(0);

    UniquePtr<FEBase> fe(FEBase::build(dim, fe_type_half));


    // Extract the FE system and DOF map, and setup the FE object.
  //~ AutoPtr<QBase> qrule = QBase::build(QGAUSS, dim, SEVENTH);
    fe->attach_quadrature_rule(qrule.get());
    const std::vector<double>& JxW = fe->get_JxW();
    const std::vector<std::vector<double> >& phi = fe->get_phi();

    // Zero out the F_g force.
    F_g.zero();

    // Loop over the local elements to compute the local integrals.
    
   

    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        const Elem* const elem = *el_it;
        fe->reinit(elem);

        const unsigned int n_qp = qrule->n_points();
        for (unsigned int qp = 0; qp < n_qp; ++qp)
        {

			 for (int d = 0; d < 3; ++d)
				F_g(d) += (rho_s - rho) * grav_const[d] * JxW[qp];
            

        }
    }
    SAMRAI_MPI::sumReduction(&F_g(0), NDIM);
    
    X_half_system.solution->close();
    

	return;
} //calculateGravitationalForce



void getSkewSymmetricAngVelTensor(TensorValue<double>& Omega,
								  VectorValue<double> W)
{

	TBOX_ASSERT( NDIM == 3); // The code is currently setup only for 3D cases //
	
	Omega.zero();

	Omega(0,1) = - W(2);
	Omega(0,2) = W(1);
	Omega(1,0) = W(2);
	Omega(1,2) = - W(0);
	Omega(2,0) = - W(1);
	Omega(2,1) = W(0);
	
	
	return;
}

void Solve6DOFSystemofEquations(const double dt, 
								VectorValue<double>& V_current,              // current linear velocity of the body
								VectorValue<double>& V_new,              // linear velocity of the body
								VectorValue<double>& W_new,              // angular velocity of the body
								VectorValue<double>& x_new,             // New position of the body
								TensorValue<double>& Q_new,                 // Rotation Matrix
								double M,							// Mass of the body
								VectorValue<double> x_com,             // New position of the body
								TensorValue<double> I_w_current,  // moment of inertia tensor of the body
								TensorValue<double>& I_w_new,  // moment of inertia tensor
								TensorValue<double> I_w_0,  // initial moment of inertia tensor
								VectorValue<double> F_b,   // total external body force   
								VectorValue<double> F_s,   // total external surface force
								VectorValue<double> T   // Torque applied on the surface
								)
{
	
	const double TOL = sqrt(std::numeric_limits<double>::epsilon());
	
	

	// This time-stepping scheme is implemented from the paper by Akkerman et al., J of Applied Mechanics,2012	

	
	TensorValue<double> Q_new_iter, I_w_new_iter, Omega_current, Omega_new, Q_current;
	VectorValue<double> W_current;
	Q_new_iter.zero();
	Omega_current.zero();
	Omega_new.zero();
	I_w_new_iter.zero();
	
	V_current = V_new;
	W_current = W_new;
	Q_current = Q_new;


//	TBOX_ASSERT((x_new - x_com).size() <= TOL && (I_w_new - I_w_current).size()<=TOL); //checking to make sure the structure center of mass coincides with the 6DOF solution
	
	
	V_new = dt * (F_b + F_s) / M + V_current;	
	x_new = 0.5 * dt * ( V_new + V_current) + x_new;
		
	int iter = 0;
	
	
	while ( (Q_new_iter - Q_new).norm() > TOL || (I_w_new_iter - I_w_new).norm() > TOL )
	{
		Q_new_iter = Q_new;
		I_w_new_iter = I_w_new;
		getSkewSymmetricAngVelTensor(Omega_current, W_current);
		getSkewSymmetricAngVelTensor(Omega_new, W_new);
		I_w_new = Q_new * I_w_0 * Q_new.transpose();
		W_new   = I_w_new.inverse() * ( dt * T + I_w_current * W_current );
		Q_new = Q_current + 0.25 * dt * (Omega_new + Omega_current) * (Q_new + Q_current);
		++iter;
	}
	
	
	pout <<"\n\n"<< "Number of 6DOF iterations = " << iter<<"\n\n";
	
	pout <<"\n\n"<< " Translationl Velocity = " <<  V_new(1) <<"\n\n";

	
    
	return;
}//Solve6DOFSystemofEquations


void updateVelocityAndPositionOfSolidPoints(VectorValue<double> x_com,
										  VectorValue<double> V_new,              // linear velocity of the body
										  VectorValue<double> V_current,              // linear velocity of the body
										  VectorValue<double> W,              // angular velocity of the body
										  const double dt,
										  EquationSystems& solid_equation_systems)
{

	
                VectorValue<double>  RR, WxR, X_new;
                
				RR.zero();
				WxR.zero();
				X_new.zero();
				
				
                MeshBase& mesh = solid_equation_systems.get_mesh();
                System& X_system = solid_equation_systems.get_system("position_new");
                const unsigned int X_sys_num = X_system.number();
                
                NumericVector<double>& X_coords = *X_system.solution;
                System& U_system = solid_equation_systems.get_system("velocity_new");
                const unsigned int U_sys_num = U_system.number();
                NumericVector<double>& U_coords = *U_system.solution;
                

                
                System& X_half_system = solid_equation_systems.get_system("position_half");
                NumericVector<double>& X_half_coords = *X_half_system.solution;
                
				

                for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
                {
                    Node* n = *it;
                    if (n->n_vars(X_sys_num))
                    {
                        TBOX_ASSERT(n->n_vars(X_sys_num) == NDIM);
                        const libMesh::Point& X = *n;
                        RR = X ;
                        WxR = W.cross(RR);
                        X_new = X + x_com + 0.5 * dt * (2.0 * WxR + (V_new + V_current));

                        
                      
                        for (unsigned int d = 0; d < NDIM; ++d)
                        {
                            const int u_dof_index = n->dof_number(U_sys_num, d, 0);
                            const int x_dof_index = n->dof_number(X_sys_num, d, 0);
                            X_coords.set(x_dof_index, X_new(d));
							U_coords.set(u_dof_index, V_new(d) + WxR(d));

                            X_half_coords.set(x_dof_index, 0.5 * (X(d) + X_new(d)));

                        }
                    }
                }
                X_coords.close();
                X_system.get_dof_map().enforce_constraints_exactly(X_system, &X_coords);
                X_system.solution->localize(*X_system.current_local_solution);
                

                
                X_half_coords.close();
                X_half_system.get_dof_map().enforce_constraints_exactly(X_half_system, &X_coords);
                X_half_system.solution->localize(*X_half_system.current_local_solution);
                
                U_coords.close();
                U_system.get_dof_map().enforce_constraints_exactly(U_system, &U_coords);
                U_system.solution->localize(*U_system.current_local_solution);
     

    return;
	
	
} //updateVelocityAndPositionOfSolidPoints


int main(int argc, char* argv[])
{
    // Initialize libMesh, PETSc, MPI, and SAMRAI.
    LibMeshInit init(argc, argv);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    


    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

 
        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();
        const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
        const string exodus_solid_filename = "solid_output.ex2"; // app_initializer->getExodusIIFilename();
        const string exodus_bndry_filename = app_initializer->getExodusIIFilename();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }


        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create a simple FE mesh.
        Mesh solid_mesh(init.comm(), NDIM);
        const double dx = input_db->getDouble("DX");
        const double ds = input_db->getDouble("MFAC") * dx;
        string elem_type = input_db->getString("ELEM_TYPE");
        const double rho = input_db->getDouble("RHO"); 
        const double rho_s = input_db->getDouble("RHO_S");

        //~ const double grav_const =input_db->getDouble("RHO");

        const double R = input_db->getDouble("R");
        if (NDIM == 2 && (elem_type == "TRI3" || elem_type == "TRI6"))
        {
#ifdef LIBMESH_HAVE_TRIANGLE
            const int num_circum_nodes = ceil(2.0 * M_PI * R / ds);
            for (int k = 0; k < num_circum_nodes; ++k)
            {
                const double theta = 2.0 * M_PI * static_cast<double>(k) / static_cast<double>(num_circum_nodes);
                solid_mesh.add_point(libMesh::Point(R * cos(theta), R * sin(theta)));
            }
            TriangleInterface triangle(solid_mesh);
            triangle.triangulation_type() = TriangleInterface::GENERATE_CONVEX_HULL;
            triangle.elem_type() = Utility::string_to_enum<ElemType>(elem_type);
            triangle.desired_area() = 1.5 * sqrt(3.0) / 4.0 * ds * ds;
            triangle.insert_extra_points() = true;
            triangle.smooth_after_generating() = true;
            triangle.triangulate();
#else
            TBOX_ERROR("ERROR: libMesh appears to have been configured without support for Triangle,\n"
                       << "       but Triangle is required for TRI3 or TRI6 elements.\n");
#endif
        }
        else
        {
            // NOTE: number of segments along boundary is 4*2^r.
            const double num_circum_segments = 2.0 * M_PI * R / ds;
            const int r = log2(0.25 * num_circum_segments);
            MeshTools::Generation::build_sphere(solid_mesh, R, r, Utility::string_to_enum<ElemType>(elem_type));
        }

        // Ensure nodes on the surface are on the analytic boundary.
        MeshBase::element_iterator el_end = solid_mesh.elements_end();
        for (MeshBase::element_iterator el = solid_mesh.elements_begin(); el != el_end; ++el)
        {
            Elem* const elem = *el;
            for (unsigned int side = 0; side < elem->n_sides(); ++side)
            {
                const bool at_mesh_bdry = !elem->neighbor(side);
                if (!at_mesh_bdry) continue;
                for (unsigned int k = 0; k < elem->n_nodes(); ++k)
                {
                    if (!elem->is_node_on_side(k, side)) continue;
                    Node& n = *elem->get_node(k);
                    n = R * n.unit();
                }
            }
        }
        solid_mesh.prepare_for_use();

        BoundaryMesh bndry_mesh(solid_mesh.comm(), solid_mesh.mesh_dimension() - 1);
        solid_mesh.boundary_info->sync(bndry_mesh);
        bndry_mesh.prepare_for_use();

        kappa_s = input_db->getDouble("KAPPA_S");
        eta_s = input_db->getDouble("ETA_S");
  
       
        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator;
        const string solver_type = app_initializer->getComponentDatabase("Main")->getString("solver_type");
        if (solver_type == "STAGGERED")
        {
            navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                "INSStaggeredHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
        }
        else if (solver_type == "COLLOCATED")
        {
            navier_stokes_integrator = new INSCollocatedHierarchyIntegrator(
                "INSCollocatedHierarchyIntegrator",
                app_initializer->getComponentDatabase("INSCollocatedHierarchyIntegrator"));
        }
        else
        {
            TBOX_ERROR("Unsupported solver type: " << solver_type << "\n"
                                                   << "Valid options are: COLLOCATED, STAGGERED");
        }
        Pointer<IBFEMethod> ib_method_ops =
            new IBFEMethod("IBFEMethod",
                           app_initializer->getComponentDatabase("IBFEMethod"),
                           &bndry_mesh,
                           app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"));
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);
        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);
        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IBFE solver.
        ib_method_ops->initializeFEEquationSystems();
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(IBFEMethod::VELOCITY_SYSTEM_NAME, vars));
		IBFEMethod::LagForceFcnData body_fcn_data(tether_force_function, sys_data);
		ib_method_ops->registerLagForceFunction(body_fcn_data);
        EquationSystems* bndry_equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        // Setup solid systems.
        libMesh::EquationSystems* solid_equation_systems(new EquationSystems(solid_mesh));
        x_new_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("position_new");
        u_new_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("velocity_new");
        x_half_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("position_half");
        u_half_solid_system = &solid_equation_systems->add_system<ExplicitSystem>("velocity_half");
        
        Order order = FIRST;
        FEFamily family = LAGRANGE;
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_new_" << d;
            x_new_solid_system->add_variable(os.str(), order, family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_new_" << d;
            u_new_solid_system->add_variable(os.str(), order, family);
        }
        
        
       for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "X_half_" << d;
            x_half_solid_system->add_variable(os.str(), order, family);
        }
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream os;
            os << "U_half_" << d;
            u_half_solid_system->add_variable(os.str(), order, family);
        }
        

        
        solid_equation_systems->init();

        // This is a horrible hack to set up the position vector.
        {
            MeshBase& mesh = solid_equation_systems->get_mesh();
            System& X_new_system = solid_equation_systems->get_system("position_new");
            const unsigned int X_new_sys_num = X_new_system.number();
            NumericVector<double>& X_new_coords = *X_new_system.solution;
            for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
            {
                Node* n = *it;
                if (n->n_vars(X_new_sys_num))
                {
                    TBOX_ASSERT(n->n_vars(X_new_sys_num) == NDIM);
                    const libMesh::Point& X = *n;
                    libMesh::Point x = X;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int dof_index = n->dof_number(X_new_sys_num, d, 0);
                        X_new_coords.set(dof_index, x(d));
                    }
                }
            }
            X_new_coords.close();
            X_new_system.get_dof_map().enforce_constraints_exactly(X_new_system, &X_new_coords);
            X_new_system.solution->localize(*X_new_system.current_local_solution);
            
            
            
            
            System& X_half_system = solid_equation_systems->get_system("position_half");
            const unsigned int X_half_sys_num = X_half_system.number();
            NumericVector<double>& X_half_coords = *X_half_system.solution;
            
            
            for (MeshBase::node_iterator it = mesh.local_nodes_begin(); it != mesh.local_nodes_end(); ++it)
            {
                Node* n = *it;
                if (n->n_vars(X_half_sys_num))
                {
                    TBOX_ASSERT(n->n_vars(X_half_sys_num) == NDIM);
                    const libMesh::Point& X = *n;
                    libMesh::Point x = X;
                    for (unsigned int d = 0; d < NDIM; ++d)
                    {
                        const int dof_index = n->dof_number(X_half_sys_num, d, 0);
                        X_half_coords.set(dof_index, x(d));
                    }
                }
            }
            X_half_coords.close();
            X_half_system.get_dof_map().enforce_constraints_exactly(X_half_system, &X_half_coords);
            X_half_system.solution->localize(*X_half_system.current_local_solution);
            
        }


        x_new_solid_system->assemble_before_solve = false;
        x_new_solid_system->assemble();

        u_new_solid_system->assemble_before_solve = false;
        u_new_solid_system->assemble();
        
        
        x_half_solid_system->assemble_before_solve = false;
        x_half_solid_system->assemble();

        u_half_solid_system->assemble_before_solve = false;
        u_half_solid_system->assemble();
        
        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        
        
        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        libMesh::UniquePtr<ExodusII_IO> exodus_solid_io(uses_exodus ? new ExodusII_IO(solid_mesh) : NULL);
        libMesh::UniquePtr<ExodusII_IO> exodus_bndry_io(uses_exodus ? new ExodusII_IO(bndry_mesh) : NULL);

        // Initialize hierarchy configuration and data on all patches.
        ib_method_ops->initializeFEData();
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        // Deallocate initialization objects.
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        
        TensorValue<double> I_w_new, I_w_0, I_w_current;
      
        VectorValue<double>  W_current, V_new,V_current, W_new, F_b, F_s, Torque;
        VectorValue<double> x_com_new, x_com;
        double M_new;


		TensorValue<double> Q_new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
		TensorValue<double> Q_current(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
        I_w_new.zero();
        I_w_current.zero();
        x_com_new.zero();
        W_new.zero();
        V_new.zero();
        V_current.zero();
        W_current.zero();
        Torque.zero();
        F_s.zero();
        F_b.zero();

        
        //****************************** Initialize RBD parameters **************************************//
       
       

	calculateGeomQuantitiesOfStructure( M_new, I_w_current, x_com, rho_s, solid_equation_systems);
	
	I_w_0 = I_w_current;
	I_w_new = I_w_current;
	x_com_new = x_com;
		
	calculateGravitationalForce(F_b, rho, rho_s, solid_equation_systems);
	

	calculateFluidForceAndTorque(F_s, Torque, x_com_new, bndry_mesh, bndry_equation_systems);

	//******************************************************************************//

     
        
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                exodus_solid_io->write_timestep(
                    exodus_solid_filename, *solid_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                exodus_bndry_io->write_timestep(
                    exodus_bndry_filename, *bndry_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
            }
        }
        
        if (SAMRAI_MPI::getRank() == 0)
        {
            
            w_new_stream.open("w_3D_new.curve", ios_base::out | ios_base::trunc);
            x_com_stream.open("x_3D_new.curve", ios_base::out | ios_base::trunc);
            v_new_stream.open("v_3D_new.curve", ios_base::out | ios_base::trunc);
            w_new_stream.precision(10);
            x_com_stream.precision(10);
            v_new_stream.precision(10);
            
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        
  

        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
			
            dt = time_integrator->getMaximumTimeStepSize();
           

           time_integrator->advanceHierarchy(dt);
            
				//****************************** RBD code **************************************//
			
			calculateGeomQuantitiesOfStructure(M_new, I_w_current, x_com, rho_s, solid_equation_systems);

			calculateGravitationalForce(F_b, rho,  rho_s, solid_equation_systems);
			calculateFluidForceAndTorque(F_s, Torque, x_com_new, bndry_mesh, bndry_equation_systems);


			Solve6DOFSystemofEquations(dt, V_new, V_current, W_new, x_com_new, Q_new, M_new,  x_com, I_w_current, I_w_new, I_w_0, F_b, F_s,Torque);
											
					
			loop_time += dt;		
		    
		    					   
			updateVelocityAndPositionOfSolidPoints(x_com, V_new, V_current, W_new, dt, *solid_equation_systems);
			
			

								  

			//******************************************************************************//

            
            
            

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";


            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            
            if (dump_viz_data && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    exodus_solid_io->write_timestep(
                        exodus_solid_filename, *solid_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                    exodus_bndry_io->write_timestep(
                        exodus_bndry_filename, *bndry_equation_systems, iteration_num / viz_dump_interval + 1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            

				if (SAMRAI_MPI::getRank() == 0)
				{
					x_com_stream << loop_time << " " << x_com_new(0) << " " << x_com_new(1) << " " << x_com_new(2) << endl;
					v_new_stream << loop_time << " " << V_new(0) << " " << V_new(1) << " " << V_new(2) << endl;
					w_new_stream << loop_time << " " << W_new(0) << " " << W_new(1) << " " << W_new(2) << endl;

					
					
				}
            
            
            
            
        }
        
         x_new_solid_system->clear();
         u_new_solid_system->clear();
		 solid_equation_systems->clear();
		 bndry_equation_systems->clear(); 


        // Close the logging streams.
        if (SAMRAI_MPI::getRank() == 0)
        {
            x_com_stream.close();
            v_new_stream.close();
            w_new_stream.close();

        }
        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];
        
        

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    return true;
} // run_example


