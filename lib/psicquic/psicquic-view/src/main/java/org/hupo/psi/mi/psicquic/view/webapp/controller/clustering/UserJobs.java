package org.hupo.psi.mi.psicquic.view.webapp.controller.clustering;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.JobStatus;
import org.hupo.psi.mi.psicquic.clustering.job.dao.DaoException;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.context.annotation.Scope;
import org.springframework.stereotype.Controller;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * Jobs that have been submitted by the user.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 1.1
 */
@Controller
@Scope( "session" )
public class UserJobs implements Serializable {

    private static final Log log = LogFactory.getLog( UserJobs.class );

    private List<ClusteringJob> currentJobs;

    @Autowired
    private ClusteringContext clusteringContext;

    public UserJobs() {
        currentJobs = new ArrayList<ClusteringJob>( );
    }

//    public String clearJobs() {
//        currentJobs.clear();
//        return "interactions";
//    }

    public List<ClusteringJob> getRefreshedJobs() {
        List<ClusteringJob> refreshedJobs = new ArrayList<ClusteringJob>( currentJobs.size() );
        for ( ClusteringJob currentJob : currentJobs ) {
            ClusteringJob job = null;
            try {
                job = clusteringContext.getDaoFactory().getJobDao().getJob( currentJob.getJobId() );
            } catch ( DaoException e ) {
                log.error( "Failed to fetch a job by jobId: " + currentJob.getJobId(), e );
            }
            refreshedJobs.add( job );
        }
        currentJobs = refreshedJobs;
        return refreshedJobs;
    }

    public List<ClusteringJob> getCurrentJobs() {
        return currentJobs;
    }

    public void setCurrentJobs( List<ClusteringJob> currentJobs ) {
        this.currentJobs = currentJobs;
    }

    public int getJobCount() {
        return currentJobs.size();
    }

    public int getCompletedJobCount() {
        return countJobsByStatus( JobStatus.COMPLETED );
    }

    public int getQueuedJobCount() {
        return countJobsByStatus( JobStatus.QUEUED );
    }

    public int getRunningJobCount() {
        return countJobsByStatus( JobStatus.RUNNING );
    }

    public int getFailedJobCount() {
        return countJobsByStatus( JobStatus.FAILED );
    }

    private int countJobsByStatus( JobStatus status ) {
        int count = 0;
        for ( ClusteringJob job : getRefreshedJobs() ) {
            if( job.getStatus().equals( status ) ) {
                count++;
            }
        }
        return count;
    }
}
