package org.hupo.psi.mi.psicquic.clustering.job.batch;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.hupo.psi.mi.psicquic.clustering.ClusteringContext;
import org.hupo.psi.mi.psicquic.clustering.job.ClusteringJob;
import org.hupo.psi.mi.psicquic.clustering.job.JobStatus;
import org.hupo.psi.mi.psicquic.clustering.job.dao.JobDao;
import org.springframework.batch.core.StepContribution;
import org.springframework.batch.core.scope.context.ChunkContext;
import org.springframework.batch.core.step.tasklet.Tasklet;
import org.springframework.batch.repeat.RepeatStatus;
import org.springframework.beans.factory.annotation.Autowired;
import psidev.psi.mi.search.SearchResult;
import psidev.psi.mi.search.Searcher;
import psidev.psi.mi.tab.model.BinaryInteraction;
import psidev.psi.mi.tab.model.builder.MitabWriterUtils;
import psidev.psi.mi.tab.model.builder.PsimiTabVersion;
import uk.ac.ebi.enfin.mi.cluster.ClusterContext;
import uk.ac.ebi.enfin.mi.cluster.Encore2Binary;
import uk.ac.ebi.enfin.mi.cluster.EncoreInteraction;
import uk.ac.ebi.enfin.mi.cluster.InteractionCluster;
import uk.ac.ebi.enfin.mi.cluster.cache.CacheStrategy;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.*;

/**
 * A tasklet that wraps the Enfin clustering algorithm.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class PsicquicClustererTasklet implements Tasklet {

    private static final Log log = LogFactory.getLog( PsicquicClustererTasklet.class );

    public static final String NEW_LINE = System.getProperty( "line.separator" );

    @Autowired
    private ClusteringContext clusteringContext;

    private String registryUrl;


    ////////////////
    // Tasklet

    public RepeatStatus execute( StepContribution contribution, ChunkContext chunkContext ) throws Exception {

        final Map<String, Object> params = chunkContext.getStepContext().getJobParameters();

        final String jobId = ( String ) params.get( "jobId" );
        if ( jobId == null ) {
            throw new IllegalArgumentException( "You must give a non null jobId" );
        }

        final JobDao jobDao = clusteringContext.getDaoFactory().getJobDao();
        final ClusteringJob job = jobDao.getJob( jobId );
        if ( job == null ) {
            log.error( "The specified jobId cannot be found in storage: " + jobId );
            // TODO how can we notify the user if its job is lost ?
            return RepeatStatus.FINISHED;
        }

        if ( log.isInfoEnabled() ) log.info( "Processing clustering job: " + jobId );

        job.setStatus( JobStatus.RUNNING );
        jobDao.update( job );

        final String miql = ( String ) params.get( "miql" );
        if ( miql == null ) {
            throw new IllegalArgumentException( "You must give a non null miql" );
        }

        final String servicesStr = ( String ) params.get( "services" );
        if ( servicesStr == null ) {
            throw new IllegalArgumentException( "You must give a non null servicesStr" );
        }
        final List<String> services = new ArrayList<String>();
        if (servicesStr.contains("|")) {
            services.addAll( Arrays.asList( servicesStr.split( "\\|" ) ) );
        } else {
            services.add( servicesStr );
        }

        if ( log.isDebugEnabled() ) log.debug( "PsicquicClustererTasklet.execute(MIQL='" + miql + "', " + services + ")" );

        final File dataLocationFile = clusteringContext.getConfig().getDataLocationFile();
//        final File cacheStorageDirectory = new File( dataLocationFile, "clustering-cache" );
//        boolean created = cacheStorageDirectory.mkdirs();
//        if ( !created ) {
//            // TODO handle error and abort process.
//        }

        ClusterContext.getInstance().setCacheStrategy( CacheStrategy.IN_MEMORY );
//        ClusterContext.getInstance().setCacheStrategy( CacheStrategy.ON_DISK );
//        ClusterContext.getInstance().getCacheManager().setCacheStorage( cacheStorageDirectory );

        // Setup up enfin clustering
        InteractionCluster iC = new InteractionCluster();
        // set registryUrl
        iC.setRegistryUrl(this.registryUrl);

        iC.addQueryAcc( miql );
        for ( String service : services ) {
            iC.addQuerySource( service );
        }

        // Run clustering
        iC.setMappingIdDbNames( "uniprotkb,irefindex,ddbj/embl/genbank,refseq,chebi" );
        iC.runService();

        // Extract MITAB lines from cluster cache
        Map<Integer, EncoreInteraction> interactionMapping = iC.getInteractionMapping();
        Encore2Binary iConverter = new Encore2Binary( iC.getMappingIdDbNames() );

        if ( log.isDebugEnabled() ) log.debug( "-------- CLUSTERED MITAB (" + interactionMapping.size() + ") --------" );

        final File jobDirectory = new File( dataLocationFile, jobId );
        boolean created = jobDirectory.mkdirs();

        if ( !created ) {
            final String msg = "Impossible to create the job directory to cluster binary interactions [miql='" +
                    miql + "', services='" + servicesStr + "', jobId='" + jobId + "']";
            log.error( msg);
        }

        if ( log.isDebugEnabled() ) log.debug( "Using Job directory: " + jobDirectory.getAbsolutePath() );

        final String mitabFilename = jobId + ".tsv";
        final File mitabFile = new File( jobDirectory, mitabFilename );

        if ( log.isDebugEnabled() ) log.debug( "MITAB file: " + mitabFile.getAbsolutePath() );
        BufferedWriter out = new BufferedWriter(new FileWriter( mitabFile ));

        for ( Map.Entry<Integer, EncoreInteraction> entry : interactionMapping.entrySet() ) {
            final EncoreInteraction ei = entry.getValue();
            final BinaryInteraction bi = iConverter.getBinaryInteraction( ei );
            final String mitab = MitabWriterUtils.buildLine(bi, PsimiTabVersion.v2_5);

            log.trace( mitab );
            out.write( mitab );
        }
        out.flush();
        out.close();

        // Build a Lucene index from MITAB
        final String luceneDirectoryName = "lucene-index";
        final File luceneDirectory = new File( jobDirectory, luceneDirectoryName );
        if ( log.isDebugEnabled() ) log.debug( "Lucene directory: " + luceneDirectory.getAbsolutePath() );

        job.setLuceneIndexLocation( luceneDirectory.getAbsolutePath() );

        boolean hasHeader = false;

        try {
            Searcher.buildIndex( luceneDirectory.getAbsolutePath(), mitabFile.getAbsolutePath(), true, hasHeader );
            final SearchResult<BinaryInteraction> result = Searcher.search( "*", luceneDirectory.getAbsolutePath(), 0, 200 );
            if ( log.isInfoEnabled() ) log.info( "Indexed " + result.getTotalCount() + " MITAB documents." );
            job.setClusteredInteractionCount( result.getTotalCount() );
        } catch ( Exception e ) {
            final String msg = "An error occured while performing the indexing of MITAB data clustered for job [miql='" +
                               miql + "', services='" + servicesStr + "', jobId='" + jobId + "']";
            log.error( msg, e );

            job.setStatus( JobStatus.FAILED );
            job.setStatusMessage( msg );
            job.setStatusException( e );
            job.setCompleted( new Date() );
            jobDao.update( job );

            return RepeatStatus.FINISHED;
        }

        // Update job repository
        job.setStatus( JobStatus.COMPLETED );
        jobDao.update( job );

        return RepeatStatus.FINISHED;
    }

    public String getRegistryUrl() {
        return registryUrl;
    }

    public void setRegistryUrl(String registryUrl) {
        this.registryUrl = registryUrl;
    }
}
