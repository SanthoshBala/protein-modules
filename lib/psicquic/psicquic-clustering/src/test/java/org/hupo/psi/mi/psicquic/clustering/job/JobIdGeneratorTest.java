package org.hupo.psi.mi.psicquic.clustering.job;

import org.hupo.psi.mi.psicquic.clustering.Service;
import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * JobIdGenerator Tester.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class JobIdGeneratorTest {

    @Test
    public void getNextJobId() throws Exception {

        // Check that the generator generated unique jobIds

        final JobIdGenerator jig = new JobIdGenerator();

        final int iteration = 1000 * 100;
        Set<String> generatedIds = new HashSet( iteration );
        for ( int i = 0; i < iteration; i++ ) {
            final String jobId = jig.getNextJobId();
            if( generatedIds.contains(jobId )) {
                Assert.fail( "JobIdGenerator did not generate unique jobId after attempt " + i  );
            }
            generatedIds.add( jobId );
        }

        Assert.assertEquals( iteration, generatedIds.size() );
    }

    @Test
    public void getNextJobId_multipleInstances() throws Exception {

        // Check that the generator generated unique jobIds

        final int iteration = 1000 * 100;
        Set<String> generatedIds = new HashSet( iteration );
        for ( int i = 0; i < iteration; i++ ) {
            JobIdGenerator jig = new JobIdGenerator();
            final String jobId = jig.getNextJobId();
            if( generatedIds.contains(jobId )) {
                Assert.fail( "JobIdGenerator did not generate unique jobId after attempt " + i  );
            }
            generatedIds.add( jobId );
        }
        Assert.assertEquals( iteration, generatedIds.size() );
    }

    @Test
    public void getJobId() throws Exception {

        // Check that the generator generated unique jobIds

        final JobIdGenerator jig = new JobIdGenerator();
        
        ClusteringJob job;
        String jobId;

        job = new ClusteringJob( "9606", Arrays.asList( new Service( "IntAct" ) ) );
        jobId = jig.generateJobId( job );
        Assert.assertNotNull( jobId );
        Assert.assertEquals( "b170829b1cd9bf67c5c1eadc56b2a2bcbd0751e0", jobId );

        // same jobid if the miql is the same when trimmed

        job = new ClusteringJob( " 9606 ", Arrays.asList( new Service( "IntAct" ) ) );
        jobId = jig.generateJobId( job );
        Assert.assertNotNull( jobId );
        Assert.assertEquals( "b170829b1cd9bf67c5c1eadc56b2a2bcbd0751e0", jobId );

        // different jobid for different services

        job = new ClusteringJob( "9606", Arrays.asList( new Service( "MINT" ) ) );
        jobId = jig.generateJobId( job );
        Assert.assertNotNull( jobId );
        Assert.assertEquals( "947ee901544dedfb5c9d83575a6d41dec633bb9b", jobId );

        // the same jobId should be generated for the same list of services, irrespectively of their order.

        job = new ClusteringJob( "9606", Arrays.asList( new Service( "IntAct" ), new Service( "MINT" ) ) );
        jobId = jig.generateJobId( job );
        Assert.assertNotNull( jobId );
        Assert.assertEquals( "67d3ffa60d5cd1a4e751f97f4f965e0bd7c4e3a2", jobId );

        job = new ClusteringJob( "9606", Arrays.asList( new Service( "MINT" ), new Service( "IntAct" ) ) );
        jobId = jig.generateJobId( job );
        Assert.assertNotNull( jobId );
        Assert.assertEquals( "67d3ffa60d5cd1a4e751f97f4f965e0bd7c4e3a2", jobId );
    }
}
