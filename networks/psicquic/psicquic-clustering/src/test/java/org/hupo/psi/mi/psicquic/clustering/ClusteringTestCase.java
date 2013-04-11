package org.hupo.psi.mi.psicquic.clustering;

import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.dao.ClusteringServiceDaoFactory;
import org.hupo.psi.mi.psicquic.clustering.job.dao.DaoException;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;

import java.util.Collection;

/**
 * Base class for clustering test case.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */

@RunWith( SpringJUnit4ClassRunner.class )
@ContextConfiguration( locations = {"classpath*:/META-INF/psicquic-clustering.spring.xml",
                                    "classpath*:/META-INF/job-clustering-test.xml"} )
public abstract class ClusteringTestCase {

    @Autowired
    @Qualifier( "inMemory" )
    private ClusteringServiceDaoFactory daoFactory;

    @Autowired
    private ClusteringContext clusteringContext;

    protected ClusteringTestCase() {
    }

    @Before
    public void beforeEachTest() throws DaoException {
        Assert.assertNotNull( clusteringContext );
        Assert.assertNotNull( daoFactory );
        Assert.assertEquals( 0, daoFactory.getJobDao().countAll() );
    }

    @After
    public void afterEachTest() throws DaoException {
        // delete all remaining jobs
        final JobDao jobDao = daoFactory.getJobDao();
        final Collection<ClusteringJob> clusteringJobs = jobDao.getAll();
        for ( ClusteringJob job : clusteringJobs ) {
            jobDao.delete( job );
        }
    }

    public ClusteringServiceDaoFactory getDaoFactory() {
        return daoFactory;
    }

    public ClusteringContext getClusteringContext() {
        return clusteringContext;
    }
}
