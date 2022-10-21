class itemTypeButton {
  constructor() {
    this.hoverYN = false;
  }
}


function defineInventoryButtons(){
  //buttons of types listeners
  for (let i = 0; i < itemTypes.length; i++){
    let typeButton = document.getElementById(itemTypes[i]);

    typeButton.addEventListener("click",function(){
      currentPage = itemTypes[i];
      displayOnlyThisType(currentPage);
    });

    // dragging
    typeButton.addEventListener("dragover", (event) => {
      event.preventDefault();
      typeButton.style.backgroundColor = "rgb(173, 216, 230)";
    });
    typeButton.addEventListener("drop", (event) => {
      event.preventDefault();
      typeButton.hoverYN = true;
    });
    typeButton.addEventListener("mouseleave",function(){
      typeButton.hoverYN = false;
      //typeButton.style.backgroundColor = "gray";
    });
  }

  //tools buttons
  document.getElementById("condenseButton").addEventListener("click",function(){
    condenseItemsPlacement();
  });
  document.getElementById("sortButton").addEventListener("click",function(){
    sortItems();
    speakText('items were sorted')
  });
  document.getElementById("newItemButton").addEventListener("click",function(){
    constructNewItem('images/items/rand-box.png', findFreeCell());
    constructNewItemOnServer(items.length - 1);
  });
}

function checkSpeach(){
  var xmlHttp = new XMLHttpRequest();
  xmlHttp.open( "GET", serverAddress+"checkSpeech", false);
  xmlHttp.send( null );
  //console.log(xmlHttp.responseText);
  if (xmlHttp.responseText != "none"){
    speakText(xmlHttp.responseText);
  }
  var xmlHttp = new XMLHttpRequest();
  xmlHttp.open( "GET", serverAddress+"setSpeech?info="+"", false);
  xmlHttp.send( null );
}

function updateInfo() {
  var xmlHttp = new XMLHttpRequest();
  xmlHttp.open( "GET", serverAddress+"getUpdates", false);
  xmlHttp.send( null );
  var serverResponse = xmlHttp.responseText;
  //console.log(serverResponse);
  addjustCharPageHeight();
  checkSpeach();
  updateItemsInfo(serverResponse);
  displayOnlyThisType(currentPage);
  placeItemsInCells();
  setTimeout(updateInfo, 1000);
}

defineInventoryButtons();
createItemIcons();
preCreateItemContextMenu();
preCreateCharImageContextMenu();
showOnlyThisCharPage(currentCharacter);

document.addEventListener('keydown', event => {
  if (event.key === 'Enter') {
    document.execCommand('insertLineBreak');
    event.preventDefault()
  }
})

document.onload = function(){updateInfo();};

/*function mousemove(event){
  console.log("pageX: ",event.pageX, 
  "pageY: ", event.pageY, 
  "clientX: ", event.clientX, 
  "clientY:", event.clientY)
}
window.addEventListener('mousemove', mousemove);*/
setTimeout(updateInfo, 100);