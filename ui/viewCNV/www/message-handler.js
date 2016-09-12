Shiny.addCustomMessageHandler("setnote",
			      function(message) {
				  $('#note').val(message.note)
			      }
			     );
