/**
 * 
 */
package org.hupo.psi.mi.psicquic.ols.client;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * @author Erik Pfeiffenberger
 *
 */
public class OntologyTerm {

    private String id;
	private String name;
	private OLSClient ols;
	private String ontology;
	private SelfDiscoveringOntologyTree sdoTree;
	
	private Set<OntologyTerm> childList;
	
	public OntologyTerm(String id,String ontology, SelfDiscoveringOntologyTree sdoTree, OLSClient ols) {
		init(id, null, ontology, sdoTree, ols);
	//	System.out.println("ONTOLOGYTERM: constructer called with id "+id );
	}
	
	public OntologyTerm(String id,String name,String ontology, SelfDiscoveringOntologyTree sdoTree, OLSClient ols) {
		init(id, name, ontology, sdoTree, ols);
	//	System.out.println("ONTOLOGYTERM: constructer called with id "+id );
	}
	
	private void init(String id,String name,String ontology, SelfDiscoveringOntologyTree sdoTree, OLSClient ols){
		this.id = id;
		this.name = name;
		this.ols = ols;
		this.ontology = ontology;
		this.sdoTree = sdoTree;

		childList = null;
		this.sdoTree.addDiscoveredTerm(this);
	}
	
	@Override
	public boolean equals(Object obj) {
		boolean equal = false;
		if(obj instanceof OntologyTerm){
			OntologyTerm term = (OntologyTerm) obj;
			
			if(this.id.equalsIgnoreCase(term.id)){
				equal = true;
			} else {
				equal = false;
			}
		}
	
		return equal;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getName() {
		
		if (name == null){
			name = ols.getTermNameByTermId(id, ontology);
		} 
		
		return name;
	}

	public OLSClient getOls() {
		return ols;
	}

	public void setOls(OLSClient ols) {
		this.ols = ols;
	}

	public Set<OntologyTerm> getChildList() {
		
		if (childList == null){
			childList = new HashSet<OntologyTerm>();

			Map<String, String>childMap = ols.getTermChildren(id, ontology);
			
			for (Map.Entry<String, String> entry : childMap.entrySet()){
				
				childList.add(new OntologyTerm(entry.getKey(), entry.getValue(), ontology, sdoTree, ols));
			}
		} 
			
		return childList;
	}

	public void setChildList(Set<OntologyTerm> childList) {
		this.childList = childList;
	}
}
