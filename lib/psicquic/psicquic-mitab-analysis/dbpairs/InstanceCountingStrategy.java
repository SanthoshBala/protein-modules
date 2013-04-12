package playground.dbpairs;

import playground.dbpairs.selectors.*;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since TODO add POM version
 */
public class InstanceCountingStrategy implements CountingStrategy {
    public boolean isShowingLinePercentage() {
        return false;
    }

    public MapSelector getIdDatabaseSelector() {
        return new IdDatabase2CountSelector();
    }

    public MapSelector getAltIdDatabaseSelector() {
        return new AltIdDatabase2CountSelector();
    }

    public MapSelector getDatabaseSelector() {
        return new Database2CountSelector();
    }

    public MapSelector getAliasSelector() {
        return new Alias2CountSelector();
    }

    public MapSelector getDatabasePairsSelector() {
        return new DatabasePair2CountSelector();
    }
}
