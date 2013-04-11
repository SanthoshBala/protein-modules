package org.hupo.psi.mi.psicquic.clustering;

import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.JobStatus;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.hupo.psi.mi.psicquic.clustering.job.dao.impl.memory.InMemoryClusteringServiceDaoFactory;
import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import java.util.Arrays;

/**
 * A playground to warm up the framework.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class JobRunnerTest extends ClusteringTestCase {

    @Test
    @Ignore
    public void runJob() throws Exception {

        Assert.assertTrue( getDaoFactory() instanceof InMemoryClusteringServiceDaoFactory );

        final JobDao jobDao = getDaoFactory().getJobDao();

        // create a job
        ClusteringJob job1 = new ClusteringJob( "brca2", Arrays.asList( new Service( "IntAct" ) ) );
        jobDao.addJob( job1.getJobId(), job1 );


        ClusteringJob job2 = new ClusteringJob( "brca2", Arrays.asList( new Service( "MINT" ), new Service( "IntAct" ) ) );
        jobDao.addJob( job2.getJobId(), job2 );

        Assert.assertEquals(JobStatus.QUEUED, job1.getStatus());
        Assert.assertEquals(JobStatus.QUEUED, job2.getStatus());

        // run the next available job
        ClusteringJob nextJob = jobDao.getNextJobToRun();
        Assert.assertNotNull( nextJob );

        Assert.assertSame( job1, nextJob );
    }
}
