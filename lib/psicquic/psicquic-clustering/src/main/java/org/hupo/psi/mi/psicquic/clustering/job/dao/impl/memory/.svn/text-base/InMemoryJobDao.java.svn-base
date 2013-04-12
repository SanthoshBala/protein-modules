package org.hupo.psi.mi.psicquic.clustering.job.dao.impl.memory;

import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.JobStatus;
import org.hupo.psi.mi.psicquic.clustering.job.dao.DaoException;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Repository;

import java.io.FileNotFoundException;
import java.util.*;

/**
 * In memory job dao.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
@Repository
@Qualifier(value = "inMemoryJobDao")
public class InMemoryJobDao implements JobDao {

    private Map<String, ClusteringJob> jobid2job = new HashMap<String, ClusteringJob>();

    @Autowired
    private ClusteringContext clusteringContext;

    public InMemoryJobDao() {
    }

    public ClusteringContext getClusteringContext() {
        return clusteringContext;
    }

    public void setClusteringContext( ClusteringContext clusteringContext ) {
        this.clusteringContext = clusteringContext;
    }

    ///////////////////
    // JobDao

    public void addJob( String jobId, ClusteringJob job ) throws DaoException {
        if ( jobId == null ) {
            throw new IllegalArgumentException( "You must give a non null jobId" );
        }
        if ( job == null ) {
            throw new IllegalArgumentException( "You must give a non null job" );
        }
        if ( jobid2job.containsKey( jobId ) ) {
            throw new IllegalStateException( "You cannot add more than once a jobId: " + jobId );
        }

        jobid2job.put( jobId, job );
    }

    public ClusteringJob getJob( String jobId ) throws DaoException {
        return jobid2job.get( jobId );
    }

    public ClusteringJob removeJob( String jobId ) throws DaoException {
        ClusteringJob job = jobid2job.remove( jobId );
        if ( job == null ) {
            throw new IllegalArgumentException( "This job could not be found: " + jobId );
        }
        return job;
    }

    /**
     * Fetches the job with the most recent completion data being in state FAILED or COMPLETED.
     *
     * @return
     */
    public ClusteringJob getLastRanJob() throws DaoException {
        ClusteringJob lastRanJob = null;
        for ( Map.Entry<String, ClusteringJob> entry : jobid2job.entrySet() ) {
            ClusteringJob job = entry.getValue();
            final JobStatus jobStatus = job.getStatus();

            if ( jobStatus != null && ( jobStatus.equals( JobStatus.FAILED )
                                        || jobStatus.equals( JobStatus.COMPLETED ) ) ) {

                if ( lastRanJob == null ) {
                    lastRanJob = job;
                } else {
                    // take the one with the most recent completion date.
                    if ( job.getCompleted().after(lastRanJob.getCompleted()) ) {
                        lastRanJob = job;
                    }
                }
            }
        }

        return lastRanJob;
    }

    /**
     * Fetches the next job with status QUEUED and the oldest creation date.
     *
     * @return
     */
    public ClusteringJob getNextJobToRun() throws DaoException {
        ClusteringJob nextJob = null;

        List<ClusteringJob> jobs = new ArrayList<ClusteringJob>(jobid2job.values());
        removeNonQueued(jobs);
        Collections.sort(jobs, new Comparator<ClusteringJob>() {
            public int compare(ClusteringJob o1, ClusteringJob o2) {
                return o1.getCreated().compareTo(o2.getCreated());
            }
        });

        if (!jobs.isEmpty()) {
            nextJob = jobs.get(0);
        }

        return nextJob;
    }

    private void removeNonQueued(List<ClusteringJob> jobs) {
        final Iterator<ClusteringJob> iterator = jobs.iterator();

        while (iterator.hasNext()) {
            ClusteringJob clusteringJob = iterator.next();
            if (JobStatus.QUEUED != clusteringJob.getStatus()) {
                iterator.remove();
            }
        }
    }

    public int getCompletedJobCount() throws DaoException {
        return countJobsByStatus( JobStatus.COMPLETED );
    }

    public int getQueuedJobCount() throws DaoException {
        return countJobsByStatus( JobStatus.QUEUED );
    }

    public int getRunningJobCount() throws DaoException {
        return countJobsByStatus( JobStatus.RUNNING );
    }

    public int getFailedJobCount() throws DaoException {
        return countJobsByStatus( JobStatus.FAILED );
    }

    private int countJobsByStatus( JobStatus status ) throws DaoException {
        int count = 0;
        for ( ClusteringJob job : getAll() ) {
            if ( job.getStatus().equals( status ) ) {
                count++;
            }
        }
        return count;
    }

    //////////////////
    // BaseDao

    public int countAll() throws DaoException {
        return jobid2job.size();
    }

    public Collection<ClusteringJob> getAll() throws DaoException {
        return new ArrayList( jobid2job.values() );
    }

    public void save( ClusteringJob job ) throws DaoException {
        // nothing to do, all in memory
    }

    public void update( ClusteringJob job ) throws DaoException {
        // nothing to do, all in memory
    }

    public void delete( ClusteringJob job ) throws DaoException {
        for ( Map.Entry<String, ClusteringJob> entry : jobid2job.entrySet() ) {
            if ( entry.getValue() == job ) {
                jobid2job.remove( entry.getKey() );
                return;
            }
        }
    }
}
