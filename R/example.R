#' Example of simulation for lazy start
#'
#' @param verbose Logical type to show or do not show messages during execution
#' @param to_plot Logical type to plot or do not plot graphical results of a simulation
#' @param seed Numeric type to set seed for a simulation, if seed = NA (by default) then it will be skipped
#' @param work_dir Working directory for a simulation, by default \code{ work_dir = getwd() }
#' @param digits Number of digits in the numeric format of pck.env environment of tugHall
#'
#' @return List of results of simulation with default values for all the parameters
#' @export
#'
#' @examples
#'
#' # it takes a time for a simulation and then it will demonstrates results, \cr
#' # so, please, wait for a while
#' simulation_example( verbose = FALSE , to_plot = FALSE )
simulation_example  <-  function( verbose = TRUE , to_plot = TRUE,
                                  seed = NA, work_dir = getwd(), digits = 5 ){

    local_dir( new = work_dir )
    pck.env$digits  =  digits

    if ( !is.na( seed ) ) set.seed( seed = seed )

    if ( verbose ) print('This code will be executed: ')
    if ( verbose ) print( simulation_example )

    # Attach packages from import list
    packages  =  list(  actuar = 'rztpois',
                        randomcoloR = 'randomColor',
                        methods = 'new',
                        stats = c('aggregate', 'rbinom', 'rexp', 'rnorm', 'runif' ),
                        stringr = c('str_length', 'str_split', 'str_sub', 'str_trim', 'str_remove'),
                        utils = c('read.delim', 'read.table', 'write.table', 'globalVariables' ),
                        grDevices = c('dev.off', 'pdf', 'rgb'),
                        graphics = c('axis', 'legend', 'lines', 'par', 'plot', 'text', 'title' ),
                        withr  =  c('local_environment', 'local_par', 'local_dir', 'local_options' ))

    for( pck in names( packages ) ){
        library( package = pck, character.only = TRUE, include.only = packages[[ pck ]])
    }

    copy_files_to_Input( files = c( 'CCDS.current.txt', 'CF.txt',
                                   'cloneinit.txt','gene_hallmarks.txt',
                                   'gene_map.txt','parameters.txt' ) ,
                         dir = 'Input' )

    define_files_names()
    define_gene_location()
    define_parameters( read_fl = TRUE , file_name = './Input/parameters.txt' )
    define_compaction_factor( read_fl = TRUE , file_name = './Input/CF.txt' )
    if ( verbose ) print_parameters()

    n_c  =  0
    repeat{
        n_c  =  n_c + 1
        if ( verbose ) print('Start simulation, please, wait for a while ... or see Sim_monitoring.txt file in working folder ')
        smlt = model( )

        if ( file.exists( pck.env$cloneoutfile ) ) break
        if ( n_c  >=  n_repeat )           break
    }
    # clones        =  pck.env$clones       # smlt[[ 1 ]]
    # onco_clones   =  pck.env$onco_clones  # smlt[[ 2 ]]

    if ( verbose ) print('Save point mutations and CNA to the files in Output folder. ')
    write_pnt_clones( pck.env$pnt_clones, file_out = 'Output/point_mutations.txt' )
    write_pnt_clones( pck.env$cna_clones, file_out = 'Output/CNA_mutations.txt' )

    if ( verbose ) print('Get data of the last simulation from cloneout file. ')
    dtst = get_flow_data( pck.env$cloneoutfile, pck.env$genefile )
    pck.env$data_avg   =  dtst$data_avg
    pck.env$data_flow  =  dtst$data_flow
    pck.env$time_max   =  dtst$time_max
    pck.env$data_last  =  dtst$data_last
    cna_mut = dtst$cna_mut
    pnt_mut   =  dtst$pnt_mut
    pnt_mut_B =  dtst$pnt_mut_B
    # onco = dtst$onco
    # hall = dtst$hall



    if ( verbose ) print('Also get VAF data and save them into file Output/VAF_data.txt ' )
    pck.env$vf = get_VAF( pnt_mut, pck.env$data_last )
    pck.env$VAF  =  get_rho_VAF( vf = pck.env$vf, rho = c( 0.0, 0.1, 0.2, 0.5, 0.7, 0.9 ) ,
                                 file_name = './Output/VAF.txt' )

    pck.env$rdr_dysf  =  get_order_of_genes_dysfunction( pnt_mut = pnt_mut,
                                                         pck.env$data_last, cna_mut,
                                        file_name = './Output/order_genes_dysfunction.txt' )

    if ( to_plot ){
        plot_order_dysfunction( pck.env$rdr_dysf , pos = c(5,800), logscale = 'y', cex = 1. )

        plot_average_simulation_data( pck.env$data_avg, pck.env$time_max )

        # Main clones
        plot_clone_evolution( pck.env$data_flow, threshold = c(0.01, 1 ), lwd = 2.0,
                              hue = c(" ", "random", "red", "orange", "yellow",
                                      "green", "blue", "purple", "pink", "monochrome")[1],
                              luminosity = c(" ", "random", "light", "bright", "dark")[4],
                              yr = NA , add_initial = TRUE, log_scale = TRUE )

        # Minor clones but large amount of them
        readline('Next? ')

        plot_clone_evolution( pck.env$data_flow, threshold = c(0.0, 0.01), lwd = 2.0,
                              hue = c(" ", "random", "red", "orange", "yellow",
                                      "green", "blue", "purple", "pink", "monochrome")[1],
                              luminosity = c(" ", "random", "light", "bright", "dark")[4],
                              yr = NA , add_initial = FALSE, log_scale = TRUE )
    }


    return( get_tugHall.Environment() )
}
