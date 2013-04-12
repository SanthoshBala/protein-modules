package org.hupo.psi.mi.psicquic.clustering.job.dao.impl.memory;

import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.ClusteringTestCase;
import org.hupo.psi.mi.psicquic.clustering.Service;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.JobIdGenerator;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.junit.Assert;
import org.junit.Test;
import org.springframework.beans.factory.annotation.Autowired;

import java.util.Arrays;

/**
 * InMemoryJobDao Tester.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class InMemoryJobDaoTest extends ClusteringTestCase {

    @Autowired
    private ClusteringContext clusteringContext;

    @Test
    public void getNextJobToRun() throws Exception {

        Assert.assertNotNull( clusteringContext );

        final JobDao jobDao = getDaoFactory().getJobDao();

        // create a job
        ClusteringJob job1 = new ClusteringJob( "9606", Arrays.asList( new Service( "IntAct" ) ) );
        jobDao.addJob( new JobIdGenerator().generateJobId( job1 ), job1 );

        Thread.sleep( 500 );

        ClusteringJob job2 = new ClusteringJob( "mouse", Arrays.asList( new Service( "MINT" ), new Service( "IntAct" ) ) );
        jobDao.addJob( new JobIdGenerator().generateJobId( job2 ), job2 );

        // run the next available job
        ClusteringJob nextJob = jobDao.getNextJobToRun();
        Assert.assertNotNull( nextJob );

        Assert.assertSame( "Got job("+ nextJob.getMiql() +") instead of job("+job1.getMiql() +")", job1, nextJob );
    }
}
