#   Ejercicio 3
#       3A)
#          Descarga ADN y análisis
#          Elija un registro de geneBank bajalo en formato fasta y gb, levantelos con Biopython y explore sus
#          atributos (record, record.id, record.name, record.description, etc.)
#          ¿Es capaz de obtener “solo” la secuencia y convertirla en un “string” para manipularla
#       3B)
#          Descargue la misma secuencia automáticamente usando Biopython.
#          Ayuda:
#          Use
#          Entrez.efetch()
#          de Biopython para acceder a GenBank
#       
#       3C)
#          Descarga secuencia de Proteína y análisis
#          Pruébelo para alguna proteína que le interesa, investigue un poco que puede hacer con el registro
#          (use tab para ayudarse)
#       3D)
#          Ahora vamos a bajarlo directamente de la web para eso usamos ExPASy (puede ver si hay
#          otras..)
#          from Bio import ExPASy
#          from Bio import SwissProt
#          handle = ExPASy.
#          get_sprot_raw
#          ("O23729")
#          record = SwissProt.read(handle)
#          handle.close()
#          print(record.entry_name)
#       3E)
#          Combine lo aprendido en la clase previas (for / while loops) y/o la información disponible en
#          biopython para cargar un conjunto de secuencias de ADN y/o proteínas. (Puede bajarlos
#          manualmente a disco y realizar código para cargarlos y/o hacer el código que los baje
#          directamente de internet). Obtenga solo las secuencias como cadenas de caracteres y guardelas
#          en una lista tal que cada elemento de la lista sea una de las secuencias