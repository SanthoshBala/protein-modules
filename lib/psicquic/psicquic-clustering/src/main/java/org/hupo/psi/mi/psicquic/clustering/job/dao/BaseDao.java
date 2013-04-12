package org.hupo.psi.mi.psicquic.clustering.job.dao;

import java.util.Collection;

/**
 * TODO document this !
 *
 * @author Samuel Kerrien (skerrien@ebi.ac.uk)
 * @version $Id$
 * @since TODO add POM version
 */
public interface BaseDao<T> {
    int countAll() throws DaoException;
    Collection<T> getAll() throws DaoException;
    void save( T object ) throws DaoException;
    void update( T object ) throws DaoException;
    void delete( T object ) throws DaoException;
}
