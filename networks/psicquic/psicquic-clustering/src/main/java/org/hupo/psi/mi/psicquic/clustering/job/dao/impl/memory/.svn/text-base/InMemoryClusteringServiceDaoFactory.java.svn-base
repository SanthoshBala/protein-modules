package org.hupo.psi.mi.psicquic.clustering.job.dao.impl.memory;

import org.hupo.psi.mi.psicquic.clustering.job.dao.ClusteringServiceDaoFactory;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class InMemoryClusteringServiceDaoFactory implements ClusteringServiceDaoFactory {

    @Autowired
    @Qualifier( "fileBasedJobDao" )
    private JobDao jobDao;

    public InMemoryClusteringServiceDaoFactory() {
    }

    public JobDao getJobDao() {
          return jobDao;
    }
}
