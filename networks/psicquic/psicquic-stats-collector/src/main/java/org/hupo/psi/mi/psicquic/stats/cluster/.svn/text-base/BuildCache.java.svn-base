package org.hupo.psi.mi.psicquic.stats.cluster;

import uk.ac.ebi.enfin.mi.cluster.ClusterContext;
import uk.ac.ebi.enfin.mi.cluster.InteractionClusterAdv;
import uk.ac.ebi.enfin.mi.cluster.cache.CacheManager;
import uk.ac.ebi.enfin.mi.cluster.cache.CacheStrategy;

import java.io.File;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 1.1
 */
public class BuildCache {

    public static void main( String[] args ) {

        final ClusterContext context = ClusterContext.getInstance();
        context.setCacheStrategy( CacheStrategy.ON_DISK );

        final CacheManager cacheManager = context.getCacheManager();
        cacheManager.setResetCache( true );

        final File cacheLocation = new File( "all_psicquic.cache" );
        cacheManager.setCacheStorage( cacheLocation );
        System.out.println( "Cache stored in: " + cacheLocation.getAbsolutePath() );

        InteractionClusterAdv iC = new InteractionClusterAdv();

        /* Query one or more IDs */
        iC.addQueryAcc( "*" );

        /* sources to query */
        /* IMEX curated databases */
//        iC.addQuerySource("DIP");
        iC.addQuerySource( "IntAct" );
        iC.addQuerySource( "MINT" );
        iC.addQuerySource( "MPact" );
        iC.addQuerySource( "MatrixDB" );
        iC.addQuerySource( "MPIDB" );
        iC.addQuerySource( "BioGrid" );
        iC.addQuerySource( "ChEMBL" );
        iC.addQuerySource( "DIP" );
        iC.addQuerySource( "InnateDB" );
        iC.addQuerySource( "MPIDB" );
        iC.addQuerySource( "Reactome" );

        long start = System.currentTimeMillis();
        iC.runService();
        long stop = System.currentTimeMillis();
        System.out.println( "Time to build: " + ((stop-start)/1000) + "s" );

        System.out.println( "#Clusters: " + cacheManager.getInteractionCache().size() );

        cacheManager.shutdown();
    }
}
