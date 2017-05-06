window.onload =  function() {

    /* Get list of plot-img elements */
    var plt_list = document.getElementsByClassName("plot-img")

    /* Get container element for modal divs */
    var source = $("#modal-handle").html();

    /* Template of modal div */
    var template = Handlebars.compile(source)

    /* For each plot-img, add a modal div */
    for (i = 0; i < plt_list.length; i++) {

        /* Get plot header */
        plot_title = plt_list[i].parentElement.parentElement.firstElementChild.innerHTML

        /* Get target for current plot */
        tar = plt_list[i].attributes["data-target"].nodeValue.substring(1);

        /* Get plot source */
        plot_src = plt_list[i].src

        data = {
            modalId: tar,
            src: plot_src,
            title: plot_title
        }

        $("#modal-container").append(template(data))
    }


}

