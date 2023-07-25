/*var chosenMll = 0;
var rnaStrings = [];

document.addEventListener('DOMContentLoaded', function () {
    document.getElementById('add-string-btn').addEventListener('click', addInputField);
    //document.getElementById('minimal-loop').addEventListener('change', addMllValue);
    //document.getElementById('start-algorithm').addEventListener('click', startNAlgorithm)
});

function addInputField() {
    var dynamicTextContainer = document.getElementById('dynamic-text-container');
    var newInput = document.createElement('input');
    newInput.type = 'text';
    newInput.name = 'content';
    newInput.placeholder = 'Write your desired RNA string';
    dynamicTextContainer.appendChild(document.createElement('br'));
    dynamicTextContainer.appendChild(newInput);
}

function addMllValue(){
    chosenMll = document.getElementById('minimal-loop').value;
}

function startNAlgorithm(){
    console.log("valore scelto: ", chosenMll);
}
document.addEventListener('DOMContentLoaded', function () {
    document.getElementById('start-algorithm').addEventListener('click', function () {
        var errorMessage = document.getElementById('error-message').value;
        if (errorMessage && errorMessage !== 'None') {
            alert(errorMessage);
            //document.getElementById('error-message').value = 'None'; // Clear the errorMessage after displaying the alert
        }
    });
});
*/