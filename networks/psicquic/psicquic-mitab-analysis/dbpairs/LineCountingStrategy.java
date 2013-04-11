package playground.dbpairs;

import playground.dbpairs.selectors.*;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since TODO add POM version
 */
public class LineCountingStrategy implements CountingStrategy {
    public boolean isShowingLinePercentage() {
        return true;
    }

    public MapSelector getIdDatabaseSelector() {
        return new IdDatabase2LineCountSelector();
    }

    public MapSelector getAltIdDatabaseSelector() {
        return new AltIdDatabase2LineCountSelector();
    }

    public MapSelector getDatabaseSelector() {
        return new Database2LineCountSelector();
    }

    public MapSelector getAliasSelector() {
        return new Alias2LineCountSelector();
    }

    public MapSelector getDatabasePairsSelector() {
        return new DatabasePair2LineCountSelector();
    }
}
