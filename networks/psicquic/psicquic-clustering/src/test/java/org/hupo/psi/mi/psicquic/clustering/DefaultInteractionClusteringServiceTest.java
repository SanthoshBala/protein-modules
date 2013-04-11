package org.hupo.psi.mi.psicquic.clustering;

import junit.framework.Assert;
import org.hupo.psi.mi.psicquic.QueryResponse;
import org.hupo.psi.mi.psicquic.clustering.job.JobStatus;
import org.hupo.psi.mi.psicquic.clustering.job.PollResult;
import org.junit.Test;
import org.springframework.beans.factory.annotation.Autowired;

import java.util.Arrays;

/**
 * DefaultInteractionClusteringService Tester.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */

public class DefaultInteractionClusteringServiceTest extends ClusteringTestCase {

    @Autowired
    private InteractionClusteringService clusteringService;

    @Test
    public void submit_poll_query() throws Exception {

        final String jobId = clusteringService.submitJob( "brca2", Arrays.asList( new Service( "IntAct" ) ) );
        PollResult pollResult = clusteringService.poll( jobId );
        Assert.assertNotNull( pollResult.getStatus() );
        Assert.assertEquals( JobStatus.QUEUED, pollResult.getStatus() );

        final int maxTry = 20;
        int tryCount = maxTry;
        do {
            tryCount--;
            Assert.assertTrue( "Failed to complete job after "+  maxTry + " attempts. Job status is: "+pollResult.getStatus(), tryCount >= 0 );

            Thread.sleep( 1000 );
            pollResult = clusteringService.poll( jobId );
        } while ( ( JobStatus.QUEUED.equals( pollResult.getStatus() ) || JobStatus.RUNNING.equals( pollResult.getStatus() ) ) );

        Assert.assertEquals( JobStatus.COMPLETED, pollResult.getStatus() );

        final String format = InteractionClusteringService.RETURN_TYPE_MITAB25;
        final QueryResponse response = clusteringService.query( jobId, "*:*", 0, 200, format );

        Assert.assertNotNull( response );

        Assert.assertNotNull( response.getResultInfo() );
        Assert.assertEquals( 200, response.getResultInfo().getBlockSize() );
        Assert.assertEquals( 0, response.getResultInfo().getFirstResult() );
        Assert.assertNotNull( response.getResultInfo().getResultType() );
        Assert.assertEquals( InteractionClusteringService.RETURN_TYPE_MITAB25, response.getResultInfo().getResultType() );

        Assert.assertEquals( 16, response.getResultInfo().getTotalResults() );

        Assert.assertNotNull( response.getResultSet() );
        Assert.assertNull( response.getResultSet().getEntrySet() );
        final String mitab = response.getResultSet().getMitab();
        Assert.assertNotNull( mitab );

        Assert.assertEquals( 16, mitab.split( "\n" ).length );
    }
}
