package org.hupo.psi.mi.psicquic.stats;

/**
 * Basic configuration holder.
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since 0.1
 */
public class Config {

    public static final String DEFAULT_SMTP_CONFIG = "/META-INF/smtp.properties";

    private String psicquicRegistryUrl;

    private String smtpConfigFile;

    private String interactionMiqlQuery;
    private String publicationMiqlQuery;

    public Config() {
    }

    public String getPsicquicRegistryUrl() {
        return psicquicRegistryUrl;
    }

    public void setPsicquicRegistryUrl( String psicquicRegistryUrl ) {
        this.psicquicRegistryUrl = psicquicRegistryUrl;
    }

    public String getSmtpConfigFile() {
        return smtpConfigFile;
    }

    public void setSmtpConfigFile( String smtpConfigFile ) {
        this.smtpConfigFile = smtpConfigFile;
    }

    public boolean hasSmtpConfigFile() {
        return smtpConfigFile != null;
    }

    public String getInteractionMiqlQuery() {
        return interactionMiqlQuery;
    }

    public void setInteractionMiqlQuery(String interactionMiqlQuery) {
        this.interactionMiqlQuery = interactionMiqlQuery;
    }

    public String getPublicationMiqlQuery() {
        return publicationMiqlQuery;
    }

    public void setPublicationMiqlQuery(String publicationMiqlQuery) {
        this.publicationMiqlQuery = publicationMiqlQuery;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "Config" );
        sb.append( "{psicquicRegistryUrl='" ).append( psicquicRegistryUrl ).append( '\'' );
        sb.append( ", smtpConfigFile='" ).append( smtpConfigFile ).append( '\'' );
        sb.append( ", interaction miqlQuery='" ).append( interactionMiqlQuery ).append( '\'' );
        sb.append( ", publication miqlQuery='" ).append( publicationMiqlQuery ).append( '\'' );
        sb.append( '}' );
        return sb.toString();
    }
}
