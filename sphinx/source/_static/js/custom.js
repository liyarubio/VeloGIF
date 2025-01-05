document.addEventListener("DOMContentLoaded", function() { 
    var cells = document.querySelectorAll(".docutils td");
    cells.forEach(function(cell) {
        // 检查单元格内容是否仅为√
        if (cell.textContent.trim() === "√") {
            var link = cell.querySelector('a');
            if (link) {
                // 为链接添加类以应用样式
                link.classList.add('green-bold');
            } else {
                // 如果没有链接，直接添加类到单元格
                cell.innerHTML = '<span class="green-bold">√</span>';
            }
        }
    });
});


document.addEventListener("DOMContentLoaded", function() { 
    var cells = document.querySelectorAll(".custom-table td");
    cells.forEach(function(cell) {
        // 检查单元格内容是否仅为√
        if (cell.textContent.trim() === "√") {
            var link = cell.querySelector('a');
            if (link) {
                // 为链接添加类以应用样式
                link.classList.add('green-bold');
            } else {
                // 如果没有链接，直接添加类到单元格
                cell.innerHTML = '<span class="green-bold">√</span>';
            }
        }
    });
});