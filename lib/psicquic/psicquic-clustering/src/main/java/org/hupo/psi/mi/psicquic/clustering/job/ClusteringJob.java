package org.hupo.psi.mi.psicquic.clustering.job;

import org.hupo.psi.mi.psicquic.clustering.Service;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.List;

/**
 * A job definition.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class ClusteringJob implements Serializable {

    private String jobId;

    /**
     * The miql query we will use to extract data from the PSICQUIC services.
     */
    private String miql;

    /**
     * The list of PSICQUIC services to cluster data from.
     */
    private List<Service> services;

    /**
     * When the job was created.
     */
    private Date created;
    private Date completed; // whether it's failed or succeeded.

    private String luceneIndexLocation;

    private JobStatus status = JobStatus.QUEUED;

    private String statusMessage;

    private Throwable statusException;

    private int clusteredInteractionCount;

    //////////////////
    // Constructors

    public ClusteringJob( String miql, List<Service> services ) {
        if ( miql == null ) {
            throw new IllegalArgumentException( "You must give a non null miql" );
        }
        this.miql = miql;

        if ( services == null ) {
            throw new IllegalArgumentException( "You must give a non null services" );
        }
        if ( services.isEmpty() ) {
            throw new IllegalArgumentException( "You must give at least one PSICQUIC service" );
        }
        this.services = new ArrayList<Service>( services );
        Collections.sort( this.services );

        this.created = new Date();

        this.jobId = new JobIdGenerator().generateJobId( this );
    }

    ///////////////////////////
    // Getters and Setters

    public String getJobId() {
        return jobId;
    }

    private void setJobId( String jobId ) {
        this.jobId = jobId;
    }

    public String getMiql() {
        return miql;
    }

    private void setMiql( String miql ) {
        this.miql = miql;
    }

    public List<Service> getServices() {
        return services;
    }

    private void setServices( List<Service> services ) {
        this.services = services;
    }

    public Date getCreated() {
        return created;
    }

    private void setCreated( Date created ) {
        this.created = created;
    }

    public JobStatus getStatus() {
        return status;
    }

    public void setStatus( JobStatus status ) {
        if ( status == null ) {
            throw new IllegalStateException( "You must give a non null status" );
        }
        this.status = status;
    }

    public String getStatusMessage() {
        return statusMessage;
    }

    public void setStatusMessage( String statusMessage ) {
        this.statusMessage = statusMessage;
    }

    public Throwable getStatusException() {
        return statusException;
    }

    public void setStatusException( Throwable statusException ) {
        this.statusException = statusException;
    }

    public Date getCompleted() {
        return completed;
    }

    public void setCompleted( Date completed ) {
        this.completed = completed;
    }

    public String getLuceneIndexLocation() {
        return luceneIndexLocation;
    }

    public void setLuceneIndexLocation( String luceneIndexLocation ) {
        this.luceneIndexLocation = luceneIndexLocation;
    }

    public int getClusteredInteractionCount() {
        return clusteredInteractionCount;
    }

    public void setClusteredInteractionCount( int clusteredInteractionCount ) {
        this.clusteredInteractionCount = clusteredInteractionCount;
    }

    ///////////////
    // Object

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "ClusteringJob" );
        sb.append( "{miql='" ).append( miql ).append( '\'' );
        sb.append( ", services=" ).append( services );
        sb.append( ", created=" ).append( created );
        sb.append( ", status=" ).append( status );
        sb.append( ", luceneIndexLocation=" ).append( luceneIndexLocation );
        sb.append( '}' );
        return sb.toString();
    }

    @Override
    public boolean equals( Object o ) {
        if ( this == o ) return true;
        if ( !( o instanceof ClusteringJob ) ) return false;

        ClusteringJob job = ( ClusteringJob ) o;

        if ( !miql.equals( job.miql ) ) return false;
        if ( !services.equals( job.services ) ) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = miql.hashCode();
        result = 31 * result + services.hashCode();
        return result;
    }
}
