double step_Verlet_C(double * forces,
                     double * forces_prev,
                     double * velocities,
                     double * velocities_new,
                     double * y,
                     double t,
                     double h,
                     double mass,
                     int n_images,
                     int n_dofs_image,
                     double (* update_field) (double, double *)
                     );
