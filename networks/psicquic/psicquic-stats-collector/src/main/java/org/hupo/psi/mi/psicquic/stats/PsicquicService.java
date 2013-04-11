package org.hupo.psi.mi.psicquic.stats;

/**
 * PSICQUIC Service bean.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class PsicquicService {

    final private String name;
    final private String soapUrl;

    private String restUrl;
    private int interactionCount;

    public PsicquicService( String name, String soapUrl ) {
        if ( name == null ) {
            throw new IllegalArgumentException( "You must give a non null name" );
        }
        if ( soapUrl == null ) {
            throw new IllegalArgumentException( "You must give a non null soapUrl" );
        }
        this.name = name;
        this.soapUrl = soapUrl;
    }

    public String getRestUrl() {
        return restUrl;
    }

    public void setRestUrl(String restUrl) {
        this.restUrl = restUrl;
    }

    public String getName() {
        return name;
    }

    public String getSoapUrl() {
        return soapUrl;
    }

    public int getInteractionCount() {
        return interactionCount;
    }

    public void setInteractionCount( int interactionCount ) {
        if ( interactionCount < 0 ) {
            throw new IllegalArgumentException( "You must give a positive count: " + interactionCount );
        }
        this.interactionCount = interactionCount;
    }

    @Override
    public boolean equals( Object o ) {
        if ( this == o ) return true;
        if ( !( o instanceof PsicquicService ) ) return false;

        PsicquicService that = ( PsicquicService ) o;

        if ( !name.equals( that.name ) ) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return name.hashCode();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder( 256 );
        sb.append( "PsicquicService[" );
        sb.append( " Name: " ).append( name ).append( ", " );
        sb.append( "SoapUrl: " ).append( soapUrl ).append( " " );
        sb.append( ']' );

        return sb.toString();
    }
}
