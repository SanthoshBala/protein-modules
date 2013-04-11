package org.hupo.psi.mi.psicquic.clustering.job.dao;

import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.springframework.stereotype.Component;
import org.springframework.stereotype.Repository;

import java.io.FileNotFoundException;

/**
 * ClusteringJob DAO.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public interface JobDao extends BaseDao<ClusteringJob> {

    void addJob( String jobId, ClusteringJob job ) throws DaoException;

    ClusteringJob getJob( String jobId ) throws DaoException;

    ClusteringJob removeJob( String jobId ) throws DaoException;

    ClusteringJob getLastRanJob() throws DaoException;

    ClusteringJob getNextJobToRun() throws DaoException;

    int getCompletedJobCount() throws DaoException;

    int getQueuedJobCount() throws DaoException;

    int getRunningJobCount() throws DaoException;

    int getFailedJobCount() throws DaoException;
}
