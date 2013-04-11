package org.hupo.psi.mi.psicquic.clustering.job.dao.impl.memory;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.JobStatus;
import org.hupo.psi.mi.psicquic.clustering.job.dao.DaoException;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Repository;

import java.io.*;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

/**
 * In memory Dao that also saves jobs on disk. If a job cannot be found in memory, then perform a disk lookup.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
@Repository
@Qualifier(value = "fileBasedJobDao")
public class FileBasedJobDao extends InMemoryJobDao {

    private static final Log log = LogFactory.getLog( FileBasedJobDao.class );

    public FileBasedJobDao() {
    }

    private File getJobFile( String jobId ) {
        final File location = getClusteringContext().getConfig().getDataLocationFile();
        return new File( location, jobId + ".dat" );
    }

    private ClusteringJob readJobFromDisk( String jobId ) throws DaoException {
        ClusteringJob job = null;

        try {
            final File file = getJobFile( jobId );
            if( file.exists() ) {
                FileInputStream fio = new FileInputStream( file );
                ObjectInputStream ois = new ObjectInputStream(fio);
                job = (ClusteringJob) ois.readObject();
            }
        } catch ( Exception e ) {
            throw new DaoException( "Could not read job from disk", e );
        }

        return job;
    }

    private void saveJobToDisk( ClusteringJob job ) throws DaoException {
        try {
            final File file = getJobFile( job.getJobId() );
            // file.deleteOnExit();
            FileOutputStream fout = new FileOutputStream( file );
            ObjectOutputStream oos = new ObjectOutputStream(fout);
            oos.writeObject( job );
            oos.close();
        } catch ( IOException e ) {
            throw new DaoException( "Error while writing job to disk", e );
        }
    }

    //////////////////
    // InMemoryJobDao

    @Override
    public void addJob( String jobId, ClusteringJob job ) throws DaoException {
        saveJobToDisk( job );
        super.addJob( jobId, job );
    }

    @Override
    public ClusteringJob getJob( String jobId ) throws DaoException {
        ClusteringJob job = super.getJob( jobId );
        if( job == null) {
            // attempt to reload from disk
            job = readJobFromDisk( jobId );
            if( job != null ) {
                // make sure it gets stored in the backing map
                super.addJob( jobId, job );
            }
        }
        return job;
    }

    @Override
    public void save( ClusteringJob job ) throws DaoException {
        saveJobToDisk( job );
        super.save( job );
    }

    @Override
    public void update( ClusteringJob job ) throws DaoException {
        saveJobToDisk( job );
        super.update( job );
    }

    @Override
    public void delete( ClusteringJob job ) throws DaoException {
        File f = getJobFile( job.getJobId() );
        boolean success = f.delete();
        if( !success ) {
            if ( log.isWarnEnabled() ) {
                log.warn( "Could not delete job file: " + job.getJobId() );
            }
        }
        super.delete( job );
    }
}
