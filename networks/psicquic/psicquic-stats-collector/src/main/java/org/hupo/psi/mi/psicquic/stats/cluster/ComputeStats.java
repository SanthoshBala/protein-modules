package org.hupo.psi.mi.psicquic.stats.cluster;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import uk.ac.ebi.enfin.mi.cluster.ClusterContext;
import uk.ac.ebi.enfin.mi.cluster.EncoreInteraction;
import uk.ac.ebi.enfin.mi.cluster.cache.CacheManager;
import uk.ac.ebi.enfin.mi.cluster.cache.CacheStrategy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Browse a cache of clustered interactions and for each PSICQUIC service represented calculates 3 things:
 * <p/>
 * 1. how many binary interactions are only seen in that service
 * <p/>
 * 2. how many binary interaction are seen multiple times in that service
 * <p/>
 * 3. how many binary interactions are seen in this service and at least one other
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class ComputeStats {
    
    private static final Log log = LogFactory.getLog( ComputeStats.class );

    public static void main( String[] args ) throws IOException {

        final ClusterContext context = ClusterContext.getInstance();
        context.setCacheStrategy( CacheStrategy.ON_DISK );
        final CacheManager cacheManager = context.getCacheManager();
        cacheManager.setResetCache( false );

        final File cacheLocation = new File( "all_psicquic.cache" );
        ClusterContext.getInstance().getCacheManager().setCacheStorage( cacheLocation );
        log.info( "Processing cache stored in: " + cacheLocation.getAbsolutePath() );

        final Map<Integer, EncoreInteraction> cache = cacheManager.getInteractionCache();
        log.info( "cluster count: " + cache.size() );

        final Map<String, ServiceStat> service2stats = new HashMap<String, ServiceStat>();

        // clusters are stored iusing a cluster id ( 1 .. cache.size() )
        final int clusterCount = cache.size();
        for ( int i = 1; i <= clusterCount; i++ ) {

            EncoreInteraction interaction = cache.get( i );

            log.info( i + "" + clusterCount );

            final Map<String, Integer> service2count = new HashMap<String, Integer>();

            for ( Map.Entry<String, List<String>> dbEntry : interaction.getExperimentToDatabase().entrySet() ) {
                // TODO why is that a list ???
                final String interactionAc = dbEntry.getKey();
                final List<String> serviceName = dbEntry.getValue();
//                log.info( "\t" + interactionAc + " -> " + serviceName );
                for ( String name : serviceName ) {

                    if ( !service2count.containsKey( name ) ) {
                        service2count.put( name, 0 );
                    }
                    int count = service2count.get( name );
                    Integer newCount = count + 1;
                    service2count.put( name, newCount );
//                    log.info( "\t\t" + name + ":" + newCount );
                }
            } // evidences

            // update the services' stats
            for ( Map.Entry<String, Integer> serviceEntry : service2count.entrySet() ) {
                String name = serviceEntry.getKey();
                int count = serviceEntry.getValue();

                if ( !service2stats.containsKey( name ) ) {
                    service2stats.put( name, new ServiceStat( name ) );
                }
                final ServiceStat stat = service2stats.get( name );
                switch ( count ) {
                    case 0:
                        throw new IllegalStateException( "service: " + name );

                    case 1:
                        if ( service2count.size() > 1 ) {
                            stat.incrementMultipleEvidenceWithOtherCount();
                        } else {
                            stat.incrementUniqueCount();
                        }
                        break;

                    default:
                        // more than 1 evidence in this service
                        if ( service2count.size() > 1 ) {
                            stat.incrementMultipleEvidenceWithOtherCount();
                        } else {
                            stat.incrementMultipleEvidenceWithinCount();
                        }
                } // evidence count per service
            } // services' stats

            if( ( i % 5000 ) == 0 ) {
                 printStats( service2stats );
            }

        } // clusters

        log.info( "" );
        log.info( "--- Final Stats ---" );
        printStats( service2stats );

        export2csv( service2stats );

        log.info( "Collected stats about " + service2stats.size() + " services" );

        ClusterContext.getInstance().getCacheManager().shutdown();
    }

    public static final String NEW_LINE = "\n";
    
    private static void export2csv( Map<String, ServiceStat> service2stats ) throws IOException {
        StringBuilder header = new StringBuilder(256);
        StringBuilder uniq = new StringBuilder(256);
        StringBuilder within = new StringBuilder(256);
        StringBuilder others = new StringBuilder(256);

        for ( Map.Entry<String, ServiceStat> serviceStatEntry : service2stats.entrySet() ) {
            final ServiceStat stat = serviceStatEntry.getValue();
            header.append( stat.getName() ).append(',');
            uniq.append( stat.getUniqueCount() ).append(',');
            within.append( stat.getName() ).append(',');
            others.append( stat.getName() ).append(',');
        }

        BufferedWriter out = new BufferedWriter( new FileWriter( new File( "stats.csv" )) );
        out.write( header.append( NEW_LINE ).toString() );
        out.write( uniq.append( NEW_LINE ).toString() );
        out.write( within.append( NEW_LINE ).toString() );
        out.write( others.append( NEW_LINE ).toString() );
        out.close();
    }

    private static void printStats( Map<String, ServiceStat> service2stats ) {
        for ( Map.Entry<String, ServiceStat> serviceStatEntry : service2stats.entrySet() ) {
            final ServiceStat stat = serviceStatEntry.getValue();
            log.info( stat );
        }
    }
}
