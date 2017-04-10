from OneBeamTwoDofs import OneBeamTwoDofs
print('Running...')

#force_filename = '/home/bowen/VIV/PostProcessing/fushixiao/KC84_T7.5/force_deviation.dat'
force_filename = 'dat/force_modified.dat'
save_directory = 'dat'
p = OneBeamTwoDofs(
    force_filename=force_filename,
    save_directory=save_directory,
    filetypes=['png'],
    # for force_mode
    mass_ratio=1.525,
    horizontal=False,
    top_tension=500,
    tension_include_boyancy=False,
    length=5,
    diameter=0.024,
    bending_stiffness=10.5,
    fluid_density=1000,
    fluid_velocity=0.268)
p.plot_crossflow_force(
    start_time=40,
    end_time=80,
    medium_step=8,
    short_step=4,
    force_min=-0.2,
    force_max=0.2,
    force_delta=0.2,
    force_std=0.2,
    #velocity_min=-8,
    #velocity_max=8,
    curvature_min=-3e-4,
    curvature_max=3e-4,
    curvature_delta=2e-4,
    curvature_std=2e-4,
    frequency_min=0,
    frequency_max=8,
    mode_number_min=1,
    mode_number_max=5,
    modal_weight_force_delta=0.2,
    #reduced_velocity_max=8,
    #period=8,
    wavelet_analysis=False,
    contour_curvature=False,
    make_animation=False)
p.plot_inline_force(
    #start_time=48.75,
    #end_time=52.5,
    start_time=40,
    end_time=80,
    medium_step=8,
    short_step=8,
    force_min=-0.2,
    force_max=0.2,
    force_delta=0.2,
    force_std=0.2,
    #velocity_min=-5,
    #velocity_max=5,
    curvature_min=-25e-5,
    curvature_max=25e-5,
    curvature_delta=25e-5,
    curvature_std=15e-5,
    frequency_min=0,
    frequency_max=8,
    mode_number_min=1,
    mode_number_max=5,
    modal_weight_force_delta=0.2,
    wavelet_analysis=False,
    contour_curvature=False,
    make_animation=False)
print('Finished!\n')
