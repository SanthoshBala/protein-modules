package playground.dbpairs.selectors;

import playground.dbpairs.Report;

import java.util.Map;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since TODO add POM version
 */
public class Database2LineCountSelector implements MapSelector {
    public Map<String, Integer> select( Report report ) {
        return report.getDatabase2lineCount();
    }
}
