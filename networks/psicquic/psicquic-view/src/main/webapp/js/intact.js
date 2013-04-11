function submitEnter(commandId, e)
{
    var keycode;
    if (window.event)
        keycode = window.event.keyCode;
    else if (e)
        keycode = e.which;
    else
        return true;

    if (keycode == 13) {
        document.getElementById(commandId).click();
        return false;
    } else
        return true;
}

function showhide(id){
    if (document.getElementById){
        obj = document.getElementById(id);
        if (obj.style.display == "none"){
            obj.style.display = "";
        } else {
            obj.style.display = "none";
        }
    }
}

function show(id){
    if (document.getElementById){
        obj = document.getElementById(id);
        if (obj.style.display == "none"){
            obj.style.display = "block";
        }
    }
}

function hide(id){
    if (document.getElementById){
        obj = document.getElementById(id);
        if (obj.style.display == "block" || obj.style.display == ""){
            obj.style.display = "none";
        }
    }
}

///////////////////////////////
// Cytoscape Web interactions

function selectMerged() {
    document.getElementById('mergeOn').style.fontWeight='bold';
    document.getElementById('mergeOff').style.fontWeight='normal';

    merged = true;
    vis.edgesMerged( merged );
}

function unselectMerged() {
    document.getElementById('mergeOn').style.fontWeight='normal';
    document.getElementById('mergeOff').style.fontWeight='bold';

    merged = false;
    vis.edgesMerged( merged );
}

function selectForceDirectedLayout() {
    document.getElementById('forceDirectedLayout').style.fontWeight='bold';
    document.getElementById('radialLayout').style.fontWeight='normal';
    document.getElementById('circleLayout').style.fontWeight='normal';
    vis.layout('ForceDirected');
}

function selectRadialLayout() {
    document.getElementById('forceDirectedLayout').style.fontWeight='normal';
    document.getElementById('radialLayout').style.fontWeight='bold';
    document.getElementById('circleLayout').style.fontWeight='normal';
    vis.layout('Radial');
}

function selectCircleLayout() {
    document.getElementById('forceDirectedLayout').style.fontWeight='normal';
    document.getElementById('radialLayout').style.fontWeight='normal';
    document.getElementById('circleLayout').style.fontWeight='bold';
    vis.layout('Circle');
}

function graphResize() {
    var heigth = $(window).height() - 330;
    var x = document.getElementById('cytoscapeweb');
    if(x!=null){
        x.style.height = heigth + 'px';
    }
}

function psicquic_selectAll(serviceCount) {
    for (i=0; i<serviceCount; i++) {
        var checkbox = document.getElementById('serviceSel_'+i);

        if (!checkbox.disabled) {
            checkbox.checked = true;
        }
    }
}

function psicquic_selectNone(serviceCount) {
    for (i=0; i<serviceCount; i++) {
        document.getElementById('serviceSel_'+i).checked = false;
    }
}

function highlightSpecies( vis ) {
    if( document.getElementById('species').style.fontWeight == 'normal' ) {
        document.getElementById('species').style.fontWeight = 'bold';

        // Create the mapper:
        var colorMapper = {
            attrName: "species",
            entries: [
                { attrValue: "Human", value: "#ff0000" },
                { attrValue: "Rat", value: "#00ff00" },
                { attrValue: "Mouse", value: "#0000ff" }
            ]
        };

        // Set the mapper to a Visual Style;
        var style = {
            nodes: {
                color: { discreteMapper: colorMapper }
            }
        };

        // Set the new style to the Visualization:
        vis.visualStyle( style );

    } else {

        document.getElementById('species').style.fontWeight = 'normal';

        var colorMapper = {
            attrName: "species",
            entries: [{ attrValue: "Human", value: "#000000" }]
        };

        // Set the mapper to a Visual Style;
        var style = {
            nodes: {
                color: { discreteMapper: colorMapper }
            }
        };

        // Set the new style to the Visualization:
        vis.visualStyle( style );
    }
}
