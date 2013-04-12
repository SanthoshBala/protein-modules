package org.hupo.psi.mi.psicquic.registry;

/**
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id: ServiceFilter.java 149 2009-08-11 11:07:37Z brunoaranda $
 */
public interface ServiceFilter {

    boolean accept(ServiceType service);
}
