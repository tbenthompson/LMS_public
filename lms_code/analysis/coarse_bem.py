from run_bem import set_params, create_fault_mesh, create_surface_mesh,\
                    bemify, boundary_conditions, assemble, constrain,\
                    solve, evaluate_surface_disp
import lms_code.lib.rep2 as rep2

if __name__ == "__main__":
    geom = rep2.load('lms_geometry')

    # Compute the BEM solution
    d = dict()
    set_params(d)

    d['degree'] = 1
    d['skip_vertices'] = 4
    d['quad_max'] = 10
    d['quad_logr'] = 10
    d['quad_oneoverr'] = 10

    create_fault_mesh(d, geom)
    create_surface_mesh(d, geom)
    bemify(d)
    boundary_conditions(d)
    assemble(d)
    constrain(d)
    solve(d)
    evaluate_surface_disp(d)
    rep2.save("bem_coarse", d)
