package playground.dbpairs.selectors;

import playground.dbpairs.Report;

import java.util.Collection;

/**
* TODO document this !
*
* @author Samuel Kerrien (skerrien@ebi.ac.uk)
* @version $Id$
* @since TODO add POM version
*/
public interface CollectionSelector {
    Collection<String> select( Report report );
}
