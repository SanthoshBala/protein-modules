/**
 * 
 */
package org.hupo.psi.mi.psicquic.ols.client;

import java.util.*;

/**
 * @author Erik Pfeiffenberger
 *
 */
public class SelfDiscoveringOntologyTree {
		
	private Map<String, OntologyTerm> alreadyDiscovered;
	
	
	public Map<String, OntologyTerm> getAlreadyDiscovered() {
		return alreadyDiscovered;
	}

	public void setAlreadyDiscovered(Map<String, OntologyTerm> alreadyDiscovered) {
		this.alreadyDiscovered = alreadyDiscovered;
	}

	OLSClient ols;
	
	String ontology;
	

	
	public SelfDiscoveringOntologyTree(String ontology) {
		this.ols = new OLSClient();
		this.ontology = ontology;
		
		alreadyDiscovered = new HashMap<String, OntologyTerm>();
	}
	
	public SelfDiscoveringOntologyTree() {
		this.ols = new OLSClient();

		alreadyDiscovered = new HashMap<String, OntologyTerm>();
	}
	
	public void addDiscoveredTerm(OntologyTerm oTerm){
		
		if(!alreadyDiscovered.containsKey(oTerm.getId())){
			alreadyDiscovered.put(oTerm.getId(), oTerm);
		}
	}
	
	public boolean isChildTag(String parentId, String childId){
		parentId = parentId.toUpperCase();
		childId = childId.toUpperCase();
		
		OntologyTerm parentTerm = alreadyDiscovered.get(parentId);
		OntologyTerm childTerm = alreadyDiscovered.get(childId);
		
		if(parentTerm == null){
			parentTerm = new OntologyTerm(parentId,ontology, this, ols);
		}
		
		if(childTerm == null){
			childTerm = new OntologyTerm(childId,ontology, this, ols);
		} 
		
		return doIsChildTag(parentTerm, childTerm);
	}

	private boolean doIsChildTag(OntologyTerm parentTerm, OntologyTerm childTerm) {
		
		
		boolean isChildTag = false;
		
		if(parentTerm.getChildList().size() > 0){
			
			for (OntologyTerm parentChild : parentTerm.getChildList()){
				if(parentChild.getId().equalsIgnoreCase(childTerm.getId())){
					isChildTag = true;
					break;
				}
			}
			
			if(!isChildTag){
				for(OntologyTerm childOfParentTag : parentTerm.getChildList()){
					if(doIsChildTag(childOfParentTag, childTerm)){
						isChildTag = true;
						break;
					}
				}
			}
			
			
		} else {
			isChildTag = false;
		}
		
		return isChildTag;
	}
	
	public Set<String> getIdByName(String text){
		Collection<String> result= ols.getTermsByName(text, ontology).values();
		
		Set<String> returnVal = new HashSet<String>();
		
		for(String r : result){
			returnVal.add(r);
		}
		
		return returnVal;
		
	}
	
	public OntologyTerm getTerm(String id){
		id = id;
		OntologyTerm term = null;
		if(alreadyDiscovered.containsKey(id)){
			term = alreadyDiscovered.get(id);
		} else {
			term = new OntologyTerm(id, ontology, this, ols);
		}
		
		return term;
	}
	
	public String getOntology() {
		return ontology;
	}
	
	public void setOntology(String ontology) {
	
		this.ontology = ontology;
	}
	
	public void clearDiscovered(){
		this.alreadyDiscovered = new HashMap<String, OntologyTerm>();
	}

}
