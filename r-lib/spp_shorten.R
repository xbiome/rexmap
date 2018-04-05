# Shorten otus that map to multiple species with equal similarity
spp_shorten = function (x, gen_sep=', ', sp_sep='/', max_sp=2, max_char=33, max_gen_ch=15,
                        max_sp_ch=30) {
   # sp_sep     = separator for species within the same genus
   # gen_sep    = separator for genuses
   # max_sp     = max. no of species to show before we omit names and just put no. {#}
   # max_char   = max. no. of characters in the final output string
   # max_gen_ch = max. no. of character in genus. some are really long
   # max_sp_ch  = max. no of characters in species name (w/o genus)
   if (is.na(x)) return(NA)
   out = ''
   sps = strsplit(x, ',')[[1]]
   i = 1
   j = 0
   last_gen = ''  # Keep track of last genus
   out_sps = ''  # Collect all species before displaying
   # If there is too many of them, just show a number
   for (s in sps) {
      sp = strsplit(s, ' ')[[1]] # Extract genus sp[1] and species sp[2]
      # Check if Genus name is too long
      # if (nchar(sp[1]) > max_gen_ch) {
      #    genus = paste0(substr(sp[1], 1, max_gen_ch), '.')
      # } else {
      #    genus = sp[1]
      # }
      if (nchar(sp[1]) > max_gen_ch) {
         sp[1] = paste0(substr(sp[1], 1, max_gen_ch), '.')
      }
      genus = sp[1]
      # Check if species name is too long
      if (nchar(sp[2]) > max_sp_ch) {
         sp[2] = paste0(substr(sp[2], 1, max_sp_ch), '.')
      }
      
      if (sp[1] != last_gen) { # New genus, put the full name
         if (i > 1) {
            if (out_sps != '') { # Write last species output
               if (j > max_sp) {
                  out_sps = paste0('sp.{', j, '}')
               }
               out = paste0(out, ' ', out_sps)
            }
            out = paste0(out, gen_sep, genus)
         } else {
            out = genus
         }
         j = 1
         out_sps = sp[2]
      } else {
         # Another of the same genus
         out_sps = paste0(out_sps, sp_sep, sp[2])
         j = j+1
         # out = paste0(out, sp_sep, sp[2])
      }
      last_gen = genus
      i = i+1
   }
   # Write last species
   if (j > max_sp) {
      out_sps = paste0('sp.{', j, '}')
   }
   out = paste0(out, ' ', out_sps)
   
   # Now check if the string is over max_char limit even with this
   if (nchar(out) > max_char) {
      # Swap the species names with just sp.
      out = paste(gsub(' .*', ' sp', strsplit(out, gen_sep)[[1]]), collapse=gen_sep)
   }
   # Cleanup brackets
   out = gsub('\\[|\\]', '', out)
   return(out)
}
