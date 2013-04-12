package org.hupo.psi.mi.psicquic.stats.cluster;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since TODO add POM version
 */
public class ServiceStat {

    String name;
    int uniqueCount;
    int multipleEvidenceWithinCount;
    int multipleEvidenceWithOtherCount;

    public ServiceStat( String name ) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public int getUniqueCount() {
        return uniqueCount;
    }

    public void incrementUniqueCount() {
        this.uniqueCount++;
    }

    public int getMultipleEvidenceWithinCount() {
        return multipleEvidenceWithinCount;
    }

    public void incrementMultipleEvidenceWithinCount() {
        this.multipleEvidenceWithinCount++;
    }

    public int getMultipleEvidenceWithOtherCount() {
        return multipleEvidenceWithOtherCount;
    }

    public void incrementMultipleEvidenceWithOtherCount() {
        this.multipleEvidenceWithOtherCount++;
    }

    @Override
    public String toString() {
        final StringBuilder sb = new StringBuilder();
        sb.append( "ServiceStat" );
        sb.append( "{name='" ).append( name ).append( '\'' );
        sb.append( ", uniqueCount=" ).append( uniqueCount );
        sb.append( ", multipleEvidenceWithinCount=" ).append( multipleEvidenceWithinCount );
        sb.append( ", multipleEvidenceWithOtherCount=" ).append( multipleEvidenceWithOtherCount );
        sb.append( '}' );
        return sb.toString();
    }
}
