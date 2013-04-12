package org.hupo.psi.mi.psicquic.wsclient.result;

import org.hupo.psi.mi.psicquic.wsclient.PsicquicClientException;
import psidev.psi.mi.tab.PsimiTabException;
import psidev.psi.mi.tab.PsimiTabReader;
import psidev.psi.mi.tab.model.BinaryInteraction;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * results containing mitab
 *
 * @author Marine Dumousseau (marine@ebi.ac.uk)
 * @version $Id$
 * @since <pre>14/08/12</pre>
 */

public class MitabSearchResult implements Serializable {

    private List<BinaryInteraction> data;
    private Integer totalCount;
    private Integer firstResult;
    private Integer maxResults;
    private String mitab;
    private PsimiTabReader mitabReader;

    public MitabSearchResult() {
        this.mitabReader = new PsimiTabReader();
    }

    public MitabSearchResult(String data, Integer totalCount, Integer firstResult, Integer maxResults)
    {
        this.mitab = data;

        this.totalCount = totalCount;
        this.firstResult = firstResult;
        this.maxResults = maxResults;
        this.mitabReader = new PsimiTabReader();
    }

    public Integer getFirstResult() {
        return firstResult;
    }

    public List<BinaryInteraction> getData() throws PsicquicClientException {
        if (data == null){
            if (mitab == null){
                return Collections.EMPTY_LIST;
            }
            else {
                try {
                    this.data = new ArrayList<BinaryInteraction>(mitabReader.read(this.mitab));
                } catch (IOException e) {
                    throw new PsicquicClientException("Problem converting the results to BinaryInteractions", e);
                } catch (PsimiTabException e) {
                    throw new PsicquicClientException("Problem converting the results to BinaryInteractions", e);
                }
            }
        }
        return data;
    }

    public Integer getMaxResults() {
        return maxResults;
    }

    public Integer getTotalCount() {
        return totalCount;
    }

    public String getRawResults(){
        return this.mitab;
    }
}
