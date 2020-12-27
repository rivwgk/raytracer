project('raytracer', 'fortran')
ifort = find_program('ifort')
add_global_arguments(['-fno-range-check', '-qopenmp', '-O3', '-xHOST'], language : 'fortran')
subdir('src.d')
target = custom_target('render',
    output : 'render',
    input : sources,
    command : [ifort, '@INPUT@', '-o', '@OUTPUT@', '-qopenmp', '-xHOST'],
    install : true,
    install_dir : 'build')
