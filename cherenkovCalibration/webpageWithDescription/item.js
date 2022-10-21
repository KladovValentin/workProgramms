class item {
    constructor() {  // Constructor
      this.itemName = "wand on intelligence";
      this.itemType = "randomStuff";
      this.itemDescription = "huy huy huy";
      this.itemWeight = 0;
      this.itemVolume = 1;
      this.itemAmount = 0;
      this.itemCell = 0;
      this.itemId = 0;
      this.equipped = false;
      this.whoseItemIs = "";
      this.itemSecondaryType = "randomStuff";
      this.itemIconSrc = "images/items/rand-box.png";
      /*this.itemDelete = itemDelete;*/
    }
}

function stringFromBool(boolValue){
  return boolValue == true ? "true" : "false";
}
function boolFromString(stringValue){
  return stringValue == "true" ? true : false;
}

function constructItemDescription(tempElement){
    let newDescription = document.createElement("div");
    newDescription.contentEditable = true;
    newDescription.classList.add("itemDescription");  //style
    newDescription.id = "tempDesc";
    
  
    var baseElement = document.getElementById("net8x8");  //net on which item is
    var cellWidth = baseElement.getBoundingClientRect().width/8.;
  
    newDescription.style.top = tempElement.getBoundingClientRect().top + 'px';
    newDescription.style.left = tempElement.getBoundingClientRect().left + tempElement.getBoundingClientRect().width - 1 +'px';
    newDescription.style.zIndex = "1";
    
    //name section
    let headerName = document.createElement("div");
    headerName.contentEditable = true;
    headerName.classList.add("itemDescription-itemHeader");
    headerName.id = "tempDescHeader";
    headerName.innerHTML = tempElement.itemName;
    headerName.style.zIndex = "1";
    newDescription.appendChild(headerName);
  
    //description itself section
    let DescDesc = document.createElement("div");
    DescDesc.contentEditable = true;
    DescDesc.classList.add("itemDescription-itemDesc");
    DescDesc.id = "tempDescDesc";
    DescDesc.innerHTML = tempElement.itemDescription;
    DescDesc.style.zIndex = "1";
    newDescription.appendChild(DescDesc);

    //amount section
    let minusAmount = document.createElement("div");
    minusAmount.contentEditable = false;
    minusAmount.classList.add("itemDescription-itemAmountButt");
    minusAmount.id = "tempMinusAmount";
    minusAmount.innerHTML = "-";
    minusAmount.style.zIndex = "1";
    minusAmount.addEventListener('click',function(){
      if (parseInt(document.getElementById("tempAmount").innerText)>=1){
        document.getElementById("tempAmount").innerText = String(parseInt(document.getElementById("tempAmount").innerText)-1);
      }
    })

    newDescription.appendChild(minusAmount);

    let amount = document.createElement("div");
    amount.contentEditable = true;
    amount.classList.add("itemDescription-itemAmount");
    amount.id = "tempAmount";
    amount.style.fontSize = "16px";
    amount.innerHTML = String(tempElement.itemAmount);
    amount.style.zIndex = "1";
    newDescription.appendChild(amount);

    let plusAmount = document.createElement("div");
    plusAmount.contentEditable = false;
    plusAmount.classList.add("itemDescription-itemAmountButt");
    plusAmount.id = "tempPlusAmount";
    plusAmount.innerHTML = "+";
    plusAmount.style.zIndex = "1";
    newDescription.appendChild(plusAmount);
    plusAmount.addEventListener('click',function(){
      document.getElementById("tempAmount").innerText = String(parseInt(document.getElementById("tempAmount").innerText)+1);
    })

    let price = document.createElement("div");
    price.contentEditable = true;
    price.classList.add("itemDescription-itemPrice");
    //price.backgroundColor = "rgb(255, 239, 100)";
    price.id = "tempPrice";
    price.style.fontSize = "16px";
    let tPrice = parseInt(tempElement.itemVolume);
    if (tPrice/1000 >= 1 && tPrice%10 == 0){
      price.innerText = String(tPrice/1000) + "k";
    }
    else if (tPrice/1000 < 1 || tPrice%10 != 0){
      price.innerText = String(tPrice);
    }
    price.style.zIndex = "1";
    newDescription.appendChild(price);


    //type section
    let descType = document.createElement("div");
    descType.contentEditable = false;
    descType.classList.add("itemDescription-itemType");
    descType.id = "tempType";
    descType.innerHTML = tempElement.itemType;
    descType.style.zIndex = "1";
    //newDescription.appendChild(descType);
    
    //secondary type section
    let secType = document.createElement("div");
    secType.contentEditable = false;
    secType.classList.add("itemDescription-itemSecType");
    secType.id = "tempSecType";
    secType.innerHTML = tempElement.itemSecondaryType;
    secType.style.zIndex = "1";
    //newDescription.appendChild(secType);


    var text = "Name: " + tempElement.itemName + "<br>";
    text = text + "Weight: " + tempElement.itemWeight.toString() + "<br>";
    text = text + "Description: " + tempElement.itemDescription + "<br>";
    //newDescription.innerHTML = text;
    newDescription.addEventListener("mouseleave", function(){
      //var newText = newDescription.innerText;
      //let newItemParameters = newText.split("\n");
      tempElement.itemName = headerName.innerText;
      tempElement.itemDescription = DescDesc.innerHTML;
      tempElement.itemType = descType.innerText;
      tempElement.itemSecondaryType = secType.innerText;
      tempElement.itemAmount = parseInt(amount.innerText);
      if (price.innerText.includes("k")){
        tempElement.itemVolume = parseInt(parseFloat(price.innerText.replace('k', ''))*1000);
      }
      else if (!price.innerText.includes("k")){
        tempElement.itemVolume = parseInt(price.innerText);
      }
      console.log(tempElement.itemName);
      updateItemInfo(findIndexInItems(tempElement.itemId));
      removeTempDescription();
    });
    document.querySelector("body").appendChild(newDescription);

    //place it inside main window
    let windowWidth = window.innerWidth;
    let windowHeight = window.innerHeight;
    if ( windowWidth - newDescription.getBoundingClientRect().right < 4){
      newDescription.style.right = windowWidth - newDescription.getBoundingClientRect().width - 4 + "px";
    }
    if ( windowHeight - newDescription.getBoundingClientRect().bottom < 4){
      newDescription.style.top = windowHeight - newDescription.getBoundingClientRect().height - 4 + "px";
    }
}
  
function findIndexInItems(tempId){
    for (let i = 0; i < items.length; i++) {
      if (items[i].itemId == tempId){
        return i;
      }
    }
    //console.log("findIndexInItems returned -1");
    return -1;
}
  
function findFirstFreeId(){
    var freeId = 0;
    var next = 0;
    while(1){
      next = 0;
      for (let i = 0; i < items.length; i++) {
        if (items[i].itemId == freeId){
          next = 1;
        }
      }
      if (next == 0){
        return freeId;
      }
      freeId +=1;
    }
}
  
function removeFromListById(tempId){
    items.splice(findIndexInItems(tempId), 1);
}
  
  
function constructNewItemOnServer(itemIndex){
    // send info about newly created item to the server
    var infoToSend = String(items[itemIndex].itemId);
    infoToSend +="-`-"+String(items[itemIndex].itemName);
    infoToSend +="-`-"+String(items[itemIndex].itemType);
    infoToSend +="-`-"+String(items[itemIndex].itemDescription);
    infoToSend +="-`-"+String(items[itemIndex].itemWeight);
    infoToSend +="-`-"+String(items[itemIndex].itemVolume);
    infoToSend +="-`-"+String(items[itemIndex].itemAmount);
    infoToSend +="-`-"+stringFromBool(items[itemIndex].equipped);
    infoToSend +="-`-"+String(items[itemIndex].whoseItemIs);
    infoToSend +="-`-"+String(items[itemIndex].itemSecondaryType);
    infoToSend +="-`-"+String(items[itemIndex].itemIconSrc);
    var xmlHttp = new XMLHttpRequest();
    xmlHttp.open( "GET", serverAddress+"genNewItem?info="+encodeURIComponent(infoToSend), false);
    xmlHttp.send( null );
}
  
function removeTempDescription(){
  if(document.getElementById("tempDesc")){
    document.getElementById("tempDesc").remove();
  }
}

function manageItemContextMenu(event, newItem){
  event.preventDefault();
  const { clientX: mouseX, clientY: mouseY } = event;
  const contextMenu = document.getElementById("itemPrimaryContext-menu");
  const deleteItemButton = document.getElementById("deleteItemButton");
  const equipItemButton = document.getElementById("equipItemButton");
  const changeIconItemButton = document.getElementById("changeIconItemButton");
  contextMenu.enactingItemId = newItem.itemId;
  console.log(newItem.equipped);
  equipItemButton.innerText = newItem.equipped ? "unequip" : "equip";

  contextMenu.classList.add("visible");
  contextMenu.style.top = `${mouseY}px`;
  contextMenu.style.left = `${mouseX}px`;
  let windowWidth = window.innerWidth;
  let windowHeight = window.innerHeight;
  if ( windowWidth - contextMenu.getBoundingClientRect().right < 4){
    contextMenu.style.right = windowWidth - contextMenu.getBoundingClientRect().width - 4 + "px";
  }
  if ( windowHeight - contextMenu.getBoundingClientRect().bottom < 4){
    contextMenu.style.top = windowHeight - contextMenu.getBoundingClientRect().height - 4 + "px";
  }

   //_________________quiting context menu listener
  document.querySelector("body").addEventListener("click",e => {
    if (e.target.offsetParent != contextMenu) {
      contextMenu.classList.remove("visible");
    }
  }, {once: true});
}
  
function constructNewItem(image,cellNumber){
    //cell - 0 to 7
    //set basic parameters of an item
    let newItem = document.createElement("img");
    newItem.classList.add("itemDesign");
    newItem.classList.add("item");
    //newItem.style.position = 'fixed';
    var baseElement = document.getElementById("net8x8");
    var cellWidth = baseElement.getBoundingClientRect().width/8.;
    var cellX = cellNumber%8;
    var cellY = (cellNumber-cellNumber%8)/8;
    newItem.style.top = baseElement.getBoundingClientRect().top + cellY * cellWidth + 'px';
    newItem.style.left = baseElement.getBoundingClientRect().left + cellX * cellWidth +'px';
    newItem.style.width = 2.5 + '%';
    //newItem.style.height = 2.5 + '%';
    newItem.src = image;
    newItem.itemName = "new item";
    newItem.itemType = currentPage;
    newItem.itemWeight = 10;
    newItem.itemVolume = 10;
    newItem.itemAmount = 1;
    newItem.itemDescription = "no description yet";
    newItem.itemCell = cellNumber;
    newItem.itemId = findFirstFreeId();
    newItem.itemSecondaryType = "random Stuff";
    newItem.whoseItemIs = currentCharacter;
    newItem.equipped = false;
    newItem.itemIconSrc = "images/items/rand-box.png";
    newItem.id = "item" + (items.length - 1).toString();
  
    //add event Listeners 
    
    //construct description on entering
    newItem.addEventListener("mouseenter", function(){
      removeTempDescription();
      constructItemDescription(newItem);
      newItem.addEventListener("mouseout", removeTempDescription);
      newItem.addEventListener('auxclick', function(e) {
        if (e.button == 1) {
          newItem.removeEventListener("mouseout", removeTempDescription);
        }
      })
      myTimeout = setTimeout(function() {
        newItem.removeEventListener("mouseout", removeTempDescription);
      }, 3000);
    });
    //context menu
    newItem.addEventListener('contextmenu', (event) => {
      manageItemContextMenu(event,newItem);
    })

    //change type when dragging an item on a changing type button 
    newItem.addEventListener('dragstart', (event) => {
      removeTempDescription();
    });
    newItem.addEventListener('dragend', (event) => {
      //if drag final point is one of the itemType buttons - change their type
      for (let i = 0; i < itemTypes.length; i++){
        let typeButton = document.getElementById(itemTypes[i]);
        if(typeButton.hoverYN == true){
          let currentPaget = currentPage;
          currentPage = itemTypes[i];
          let sameName = false;
          for (let j = 0; j < items.length; j++){
            if (items[j].itemType == currentPage && items[j].whoseItemIs == currentCharacter && items[j].itemName == newItem.itemName){
              items[j].itemAmount += newItem.itemAmount;
              sameName = true;
              removeItem(newItem);
              updateItemInfo(findIndexInItems(items[j].itemId));
            }
          }
          if (!sameName){
            //place item in free spot on new tab because they can have the same cell number on different pages with another item
            if (!newItem.equipped){  
              newItem.itemCell = findFreeCell();
            }
            newItem.itemType = itemTypes[i];
            console.log(newItem.itemType);
            updateItemInfo(findIndexInItems(newItem.itemId));
          }
          displayOnlyThisType(currentPage);
        }
      }

      //if drag final point is one of the character equipment sockets - change their subtype
      let charpage = document.getElementById(newItem.whoseItemIs + "CharPage");
      for (let i = 0; i < charpage.socketsIds.length; i++){
        let equipmentSocket = document.getElementById(newItem.whoseItemIs + "CharPage" + charpage.socketsIds[i]);
        if(equipmentSocket.hoverYN == true){
          newItem.itemSecondaryType = charpage.socketsTypes[i];
          if (!newItem.equipped){
            equipUnequipItem(findIndexInItems(newItem.itemId));
            break;
          }
          updateItemInfo(findIndexInItems(newItem.itemId));
          break;
        }
      }
      for (let i = 0; i < 9; i++){
        let quickBarSocket = document.getElementById("quickBarCell" + i.toString());
        if(quickBarSocket.hoverYN == true){
          newItem.itemSecondaryType = "quickBar";
          if (!newItem.equipped){
            equipUnequipItem(findIndexInItems(newItem.itemId));
            break;
          }
          updateItemInfo(findIndexInItems(newItem.itemId));
          break;
        }
      }

      //if drag final point is one of the character icon buttons - change their owner
      for (let i = 0; i < users.length; i++){
        let charbutton = document.getElementById(users[i] + "CharPage" + "IconButton");
        if(charbutton.hoverYN == true){
          let sameName = false;
          for (let j = 0; j < items.length; j++){
            if (items[j].itemType == currentPage && items[j].whoseItemIs == users[i] && items[j].itemName == newItem.itemName && !items[j].equipped){
              items[j].itemAmount += newItem.itemAmount;
              sameName = true;
              removeItem(newItem);
              updateItemInfo(findIndexInItems(items[j].itemId));
            }
          }
          if (!sameName){
            currentCharacter = users[i];
            if (!newItem.equipped){  
              newItem.itemCell = findFreeCell();
            }
            currentCharacter = newItem.whoseItemIs;
            newItem.whoseItemIs = users[i];
            updateItemInfo(findIndexInItems(newItem.itemId));
            displayOnlyThisType(currentPage);
          }
        }
      }
    })
    newItem.addEventListener('click',function(){
      if (newItem.itemSecondaryType == "quickBar" && newItem.equipped){
        let infoToSend = newItem.whoseItemIs + " wants to use " + newItem.itemName;
        var xmlHttp = new XMLHttpRequest();
        xmlHttp.open( "GET", serverAddress+"setSpeech?info="+encodeURIComponent(infoToSend), false);
        xmlHttp.send( null );
      }
    });

    
    let amountText = document.createElement("div");
    amountText.innerText = String(newItem.itemAmount);
    amountText.style.fontSize = "10px";
    amountText.style.position = "absolute";
    amountText.style.bottom = "5%";
    amountText.style.right = "5%";
    amountText.style.zIndex = "1000";
    newItem.appendChild(amountText);

    //append webpage with new item and push back items
    items.push(newItem);
    var index = items.length - 1;
    document.querySelector("body").appendChild(items[index]);

    
}

function placeItemInCell(newItem){
  newItem.style.transitionDuration = "0s";
  if(newItem.equipped == true){
    let charpage = document.getElementById(newItem.whoseItemIs + "CharPage");
    for (let i = 0; i < charpage.socketsIds.length; i++){
      let boxOfEquipmentExtra = document.getElementById(newItem.whoseItemIs + "CharPage" + charpage.socketsIds[i]);
      if (boxOfEquipmentExtra.occupationItemId == newItem.itemId){
        newItem.style.top = boxOfEquipmentExtra.getBoundingClientRect().top + 'px';
        newItem.style.left = boxOfEquipmentExtra.getBoundingClientRect().left +'px';
        boxOfEquipmentExtra.style.zIndex = "0";
        //console.log("I AM HERERERERERERER ");
        //newItem.style.transitionDuration = "0.2s";
        return;
      }
    }
    let charPageName = newItem.whoseItemIs + "CharPage";
    let boxid = transformItemSecTypeToBoxId(charPageName, newItem);
    if (newItem.itemSecondaryType == "quickBar"){
      for (let i = 0; i < 9; i++){
        if (charpage.quickBarItemIds[i] == newItem.itemId){
          boxid = "quickBarCell" + i.toString();
        }
      }
    }
    let boxOfEquipment = document.getElementById(boxid);
    //console.log(newItem.style.transitionDuration);
    newItem.style.width = boxOfEquipment.getBoundingClientRect().width + 'px';
    newItem.style.top = boxOfEquipment.getBoundingClientRect().top + 'px';
    newItem.style.left = boxOfEquipment.getBoundingClientRect().left +'px';
    return;
  }
  var baseElement = document.getElementById("net8x8");
  var cellWidth = baseElement.getBoundingClientRect().width/8.;
  var cellNumber = newItem.itemCell;
  var cellX = cellNumber%8;
  var cellY = (cellNumber-cellNumber%8)/8;
  newItem.style.top = baseElement.getBoundingClientRect().top + cellY * cellWidth + 'px';
  newItem.style.left = baseElement.getBoundingClientRect().left + cellX * cellWidth +'px';
  newItem.style.width = 2.5 + '%';
  //newItem.style.width = cellWidth + 'px';
  //newItem.style.height = cellWidth + 'px';
  newItem.style.transitionDuration = "0.2s";
  return;
}
  
function placeItemsInCells(){
  for (let i = 0; i < items.length; i++) {
    //items[i].style.transitionDuration = "0s";
    placeItemInCell(items[i]);
    //items[i].style.transitionDuration = "0.2s";
  }
}
  
function sortItems(){
    var sortTags = [];
    var sortItemsIndexes = [];
    //choose which are on this page
    for (let i = 0; i < items.length; i++) {
      if (items[i].itemType == currentPage && items[i].whoseItemIs == currentCharacter && items[i].equipped == false){
        sortTags.push(items[i].itemName);
        sortItemsIndexes.push(i);
      }
    }
    //sort by name first
    for (let i = 0; i < sortTags.length; i++) {
      for (let j = 0; j < sortTags.length - i -1; j++) {
        if (sortTags[j] > sortTags[j + 1]) {
          var tp = sortTags[j];
          sortTags[j] = sortTags[j+1];
          sortTags[j+1] = tp;
          
          tp = sortItemsIndexes[j];
          sortItemsIndexes[j] = sortItemsIndexes[j+1];
          sortItemsIndexes[j+1] = tp;
        }
      }
    }
    //sort cells
    for (let i = 0; i < sortTags.length; i++) {
      for (let j = 0; j < sortTags.length - i -1; j++) {
        if (items[sortItemsIndexes[j]].itemCell > items[sortItemsIndexes[j+1]].itemCell) {
          tp = items[sortItemsIndexes[j]].itemCell;
          items[sortItemsIndexes[j]].itemCell = items[sortItemsIndexes[j+1]].itemCell;
          items[sortItemsIndexes[j+1]].itemCell = tp;
        }
      }
    }
    placeItemsInCells();
}
  
function condenseItemsPlacement(){
    var tempCell = 0;
    //first, give places to ones that are in inventory and on the proper page
    for (let i = 0; i < items.length; i++) {
      if (items[i].equipped == false && items[i].itemType == currentPage && items[i].whoseItemIs == currentCharacter){
        items[i].itemCell = tempCell;
        tempCell+=1;
      }
    }
    sortItems();
}
  
function displayOnlyThisType(displayType){
    for (let i = 0; i < itemTypes.length; i++){
      document.getElementById(itemTypes[i]).style.backgroundColor = "gray";
    }
    document.getElementById(displayType).style.backgroundColor = "lightblue";
    for (let i = 0; i < items.length; i++) {
      if ((items[i].itemType == displayType || items[i].equipped == true) && items[i].whoseItemIs == currentCharacter){
        items[i].style.display = "inline";
      }
      else{
        items[i].style.display = "none";
      }
    }
}
  
function findFreeCell(){
    var filledCells = [];
    for(let i = 0; i < 64; i++){
      filledCells.push(0);
    }
    for (let i = 0; i < items.length; i++) {
      if (items[i].itemType != currentPage || items[i].equipped == true || items[i].whoseItemIs != currentCharacter){
        continue;
      }
      filledCells[items[i].itemCell] = 1;
    }
    for (let i = 0; i < 64; i++) {
      if(filledCells[i] == 0){
        return i;
      }
    }
}
  
function updateItemInfo(itemIndex){
    //send info about item to the server after changing name, description etc
    var infoToSend = String(items[itemIndex].itemId);
    infoToSend +="_"+String(items[itemIndex].itemId);
    infoToSend +="-`-"+String(items[itemIndex].itemName);
    infoToSend +="-`-"+String(items[itemIndex].itemType);
    infoToSend +="-`-"+String(items[itemIndex].itemDescription);
    infoToSend +="-`-"+String(items[itemIndex].itemWeight);
    infoToSend +="-`-"+String(items[itemIndex].itemVolume);
    infoToSend +="-`-"+String(items[itemIndex].itemAmount);
    infoToSend +="-`-"+stringFromBool(items[itemIndex].equipped);
    infoToSend +="-`-"+String(items[itemIndex].whoseItemIs);
    infoToSend +="-`-"+String(items[itemIndex].itemSecondaryType);
    infoToSend +="-`-"+String(items[itemIndex].itemIconSrc);
    var xmlHttp = new XMLHttpRequest();
    console.log(infoToSend);
    xmlHttp.open( "GET", serverAddress+"updateItemParameters?info="+encodeURIComponent(infoToSend), false);
    xmlHttp.send( null );
}

function setItemParameters(itemIndex,parametersArray){
  if (items[itemIndex].equipped != boolFromString(parametersArray[7])){
    items[itemIndex].itemCell = findFreeCell();
  }

  items[itemIndex].itemId = parseInt(parametersArray[0]);
  items[itemIndex].itemName = String(parametersArray[1]);
  items[itemIndex].itemType = String(parametersArray[2]);
  items[itemIndex].itemDescription = String(parametersArray[3]);
  items[itemIndex].itemWeight = parseFloat(parametersArray[4]);
  items[itemIndex].itemVolume = parseFloat(parametersArray[5]);
  items[itemIndex].itemAmount = parseInt(parametersArray[6]);
  items[itemIndex].equipped = boolFromString(parametersArray[7]);
  items[itemIndex].whoseItemIs = String(parametersArray[8]);
  items[itemIndex].itemSecondaryType = String(parametersArray[9]);
  items[itemIndex].src = String(parametersArray[10]);
  items[itemIndex].itemIconSrc = String(parametersArray[10]);

}

function updateItemsInfo(newParString){
    //get updates on items from server, just update all and add if no needed id exists
    var serverIds = [];
    let infoArrayItems = newParString.split("_");
    for(let i = 0; i < infoArrayItems.length; i++){
      let infoArrayParameters = infoArrayItems[i].split("-`-");
      var serverItemIndex = findIndexInItems(parseInt(infoArrayParameters[0]));
  
      //update if exists
      if(serverItemIndex != -1 && infoArrayParameters.length == 11){
        setItemParameters(serverItemIndex,infoArrayParameters);
        serverIds.push(parseInt(infoArrayParameters[0]));
      }
  
      //add if exists on server and not on page
      else if (serverItemIndex == -1 && infoArrayParameters.length == 11){
        let currentPaget = currentPage;
        currentPage = infoArrayParameters[2];
        let currentChart = currentCharacter;
        currentCharacter = infoArrayParameters[8];
        constructNewItem('images/items/rand-box.png', findFreeCell());
        currentCharacter = currentChart;
        currentPage = currentPaget;
        setItemParameters(items.length - 1,infoArrayParameters);
        items[items.length - 1].style.transitionDuration = "0.0s";
        //transfer to equipment if equipment
        if(items[items.length - 1].equipped==true){
          let charPageName = items[items.length - 1].whoseItemIs + "CharPage";
          let boxid = transformItemSecTypeToBoxId(charPageName, items[items.length - 1]);
          if (items[items.length - 1].itemSecondaryType == "quickBar"){
            let charpage = document.getElementById(charPageName);
            let foundYN = false;
            for (let i = 0; i < 9; i++){
              if (charpage.quickBarItemIds[i] == items[items.length - 1].itemId){
                boxid = "quickBarCell" + i.toString();
                foundYN = true;
              }
            }
            if (!foundYN){
              items[items.length - 1].equipped=false;
              equipUnequipItemBase(items.length - 1);
              boxid = "none";
            }
          }
          else if (boxid!="none"){
            if (document.getElementById(boxid).occupied == false){
              items[items.length - 1].equipped=false;
              equipUnequipItemBase(items.length - 1);
            }
          }
        }
        items[items.length - 1].style.transitionDuration = "0.2s";
        serverIds.push(parseInt(infoArrayParameters[0]));
      }
    }
  
    //remove from page if exists only at page
    for (let i = 0; i < items.length; i++) {
      var removeYN = 1;
      for (let j = 0; j < serverIds.length; j++) {
        if (serverIds[j] == items[i].itemId){
          removeYN = 0;
        }
      }
      if (removeYN == 1){
        removeItem(items[i]);
      }
    }
}