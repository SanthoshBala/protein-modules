package org.hupo.psi.mi.psicquic.clustering.job;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class PollResult {

    private JobStatus status;

    private String message;

    public PollResult( JobStatus status ) {
        if ( status == null ) {
            throw new IllegalArgumentException( "You must give a non null status" );
        }
        this.status = status;
    }

    public JobStatus getStatus() {
        return status;
    }

    public void setStatus( JobStatus status ) {
        this.status = status;
    }

    public String getMessage() {
        return message;
    }

    public void setMessage( String message ) {
        this.message = message;
    }

    /////////////////
    // Object

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder( 512 );
        sb.append( "PollResult" );
        sb.append( "{status=" ).append( status );
        sb.append( ", message='" ).append( message ).append( '\'' );
        sb.append( '}' );
        return sb.toString();
    }

    @Override
    public boolean equals( Object o ) {
        if ( this == o ) return true;
        if ( !( o instanceof PollResult ) ) return false;

        PollResult that = ( PollResult ) o;

        if ( status != that.status ) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return status.hashCode();
    }
}
