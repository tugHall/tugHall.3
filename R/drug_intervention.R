


#'  Function to emulate drug intervention to pool of tumor cells by killing a part of these cells
#'  and blocking a part of clones corresponding to malfunctioned gene
#'
#' @param kill_prob Probability of killing cancer cells corresponding to the malfunctioned gene
#' @param block_prob Probability of blocking cancer cells corresponding to the malfunctioned gene
#' @param gene Name of target gene to kill and block tumor cells by a drug
#'
#' @return  NULL changing clones and onco_clones objects in tugHall environment pck.env
#' @export
#'
#' @examples
drug_intervention  <-  function( kill_prob = 0, block_prob = 1, gene ){

    local_environment( env = pck.env )

    n_gene  =  which( pck.env$onco$name == gene )
    if ( length( n_gene ) == 0 ) stop( 'There is no defined gene in the pool of clones with genes of interest')

    slct  =  which ( sapply( X = 1:length( pck.env$clones ),
                        FUN = function( x ) pck.env$clones[[ x ]]$gene[ n_gene ] ==  1 ) )
    if (length( slct ) > 0 ){
        rest  =   ( 1 : length( pck.env$clones ))[ - slct ]

        # Kill the cells by calc_binom( 1, N_cells, (1 - kill_prob) ) function:
        sapply( X = slct, FUN = function( x) pck.env$clones[[ x ]]$N_cells  =
                    calc_binom( 1, pck.env$clones[[ x ]]$N_cells, (1 - kill_prob) ) )

        # Save copies of clones (splitting the blocked clones):
        clones_to_change  =   sapply( X = slct, FUN = function( x ) clone_copy( pck.env$clones[[ x ]] ) )
        onco_clones_to_change  =   sapply( X = slct, FUN = function( x ) onco_copy( pck.env$onco_clones[[ x ]] ) )
        sapply( X = 1:length( slct ), FUN = function( x ) onco_clones_to_change[[ x ]]$id  =  clones_to_change[[ x ]]$id )

        # Number of blocking cells for each clone:

        N_block  =  sapply( X = slct, FUN = function( x) calc_binom( 1, pck.env$clones[[ x ]]$N_cells , block_prob ) )

        sapply( X = 1:length( slct ), FUN = function( x ) {
                    pck.env$clones[[ slct[x] ]]$N_cells  = pck.env$clones[[ slct[x] ]]$N_cells - N_block[ x ]
        } )
        sapply( X = 1:length( slct ), FUN = function( x ) clones_to_change[[ x ]]$N_cells  =  N_block[ x ] )

    }



    # pck.env$clones  =  sapply(X = rest, function( x ) pck.env$clones[[ x ]] )

    # pck.env$onco_clones  =  sapply(X = rest, function( x ) pck.env$onco_clones[[ x ]] )


    return( NULL )
}

