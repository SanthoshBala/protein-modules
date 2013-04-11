package org.hupo.psi.mi.psicquic.clustering;

import org.junit.Assert;
import org.junit.Test;
import org.springframework.batch.core.launch.JobLauncher;
import org.springframework.beans.factory.annotation.Autowired;

/**
 * SpringDependencies Tester.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class SpringDependenciesTest extends ClusteringTestCase {

    @Autowired
    private ClusteringContext clusteringContext;

    @Autowired
    private InteractionClusteringService clusteringService;

    @Autowired
    private ClusteringLauncher launcher;

    @Test
    public void context() throws Exception {
        Assert.assertNotNull( clusteringContext );

        Assert.assertNotNull( clusteringContext.getConfig() );
        Assert.assertNotNull( clusteringContext.getDaoFactory() );
        Assert.assertNotNull( clusteringContext.getSpringContext() );
    }

    @Test
    public void clusteringService() throws Exception {
        Assert.assertNotNull( clusteringService );
        Assert.assertNotNull( clusteringService.getClusteringContext() );
    }

    @Test
    public void launcher() throws Exception {
        Assert.assertNotNull( launcher );

        Assert.assertNotNull( launcher.getClusteringContext() );
        Assert.assertNotNull( launcher.getJob() );
        Assert.assertNotNull( launcher.getJobLauncher() );
        Assert.assertNotNull( launcher.getJobRepository() );
    }
}
