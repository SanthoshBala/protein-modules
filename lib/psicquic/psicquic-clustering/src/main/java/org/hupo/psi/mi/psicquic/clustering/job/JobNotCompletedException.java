package org.hupo.psi.mi.psicquic.clustering.job;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class JobNotCompletedException extends Exception {
    public JobNotCompletedException() {
        super();
    }

    public JobNotCompletedException( String message ) {
        super( message );
    }

    public JobNotCompletedException( String message, Throwable cause ) {
        super( message, cause );
    }

    public JobNotCompletedException( Throwable cause ) {
        super( cause );
    }
}
