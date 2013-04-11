package org.hupo.psi.mi.psicquic.clustering;

import org.hupo.psi.mi.psicquic.QueryResponse;
import org.hupo.psi.mi.psicquic.clustering.job.PollResult;

import java.util.Arrays;
import java.util.List;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public interface InteractionClusteringService {

    public static final int BLOCKSIZE_MAX = 200;

    public static final String RETURN_TYPE_XML25 = "psi-mi/xml25";
    public static final String RETURN_TYPE_MITAB25 = "psi-mi/tab25";
    public static final String RETURN_TYPE_COUNT = "count";

    public static final String RETURN_TYPE_DEFAULT = RETURN_TYPE_MITAB25;

    public static final List<String> SUPPORTED_RETURN_TYPES = Arrays.asList( RETURN_TYPE_XML25, RETURN_TYPE_MITAB25, RETURN_TYPE_COUNT );

    ClusteringContext getClusteringContext();

    void setClusteringContext( ClusteringContext clusteringContext );

    /**
     * Initiate the clustering of a MIQL query based on the given list of services.
     * <br/>
     * This job is asynchronous and the user should call <code>poll( jobId )</code> to find out the status of the job.
     *
     * @param miql
     * @param services
     * @return the jobId assigned.
     */
    String submitJob( String miql, List<Service> services );

    /**
     * Provides a status to a given job, identified by the given jobId.
     *
     * @param jobId
     * @return a poll result that gives indication of the status of the job and optionaly a message.
     */
    PollResult poll( String jobId );

    /**
     * Once a job has a status <code>JobStatus.COMPLETED</code>, a user can query the clustered data.
     * @param jobId
     * @param query run a query on the clustered data, * will return all.
     * @param from considering the data being indexed from 0..(n-1), we return data from the interaction position at index 'from'
     * @param maxResult the maximum number of interactions to be returned.
     * @param resultType the return type of the data generated.
     * @return  a <code>QueryResponse</code> wrapping the resulting data.
     */
    QueryResponse query( String jobId,
                         String query,
                         final int from,
                         final int maxResult,
                         String resultType ) throws ClusteringServiceException;
}
