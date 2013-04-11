package org.hupo.psi.mi.psicquic.clustering;

import java.io.Serializable;

/**
 * Basic description of a PSICQUIC service.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class Service implements Comparable<Service>, Serializable {

    String name;
    String restUrl;

    public Service( String name ) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setName( String name ) {
        this.name = name;
    }

    public String getRestUrl() {
        return restUrl;
    }

    public void setRestUrl( String restUrl ) {
        this.restUrl = restUrl;
    }

    //////////////////
    // Comparable

    public int compareTo( Service s ) {
        return this.getName().compareTo( s.getName() );
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "Service" );
        sb.append( "{name='" ).append( name ).append( '\'' );
        sb.append( ", restUrl='" ).append( restUrl ).append( '\'' );
        sb.append( '}' );
        return sb.toString();
    }
}
