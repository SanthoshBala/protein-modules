package org.hupo.psi.mi.psicquic.view.webapp.controller;

/**
 * Search webapp exception.
 *
 * @author Bruno Aranda (baranda@ebi.ac.uk)
 * @version $Id: SearchWebappException.java 1482 2012-10-18 11:18:29Z noedelta $
 */
public class SearchWebappException extends RuntimeException{
    public SearchWebappException() {
    }

    public SearchWebappException(String s) {
        super(s);
    }

    public SearchWebappException(String s, Throwable throwable) {
        super(s, throwable);
    }

    public SearchWebappException(Throwable throwable) {
        super(throwable);
    }
}
