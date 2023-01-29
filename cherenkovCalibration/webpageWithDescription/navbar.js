
function makeLink(title, href){
    let homePage = document.createElement('a');
    homePage.appendChild(document.createTextNode(title));
    homePage.title = title;
    homePage.href =  href;
    return homePage
}   


let navbar = document.createElement("div");
navbar.classList.add("navbar"); 
navbar.id = "tempNav";

navbar.appendChild(makeLink("Home", "index.html"));