use bigtools;
use bigtools::beddata::BedParserStreamingIterator;
use std::thread::Builder;
use std::collections::HashMap;
use std::path::Path;

use bigtools::BigWigWrite;


/// Represents a single value in a bigWig file
#[derive(Copy, Clone, Debug, PartialEq)]
//#[cfg_attr(feature = "write", derive(Serialize, Deserialize))]
pub struct Value {
    pub start: u32,
    pub end: u32,
    pub value: f32,
}

impl Value{
    pub fn flat(&self) -> ( u32, u32, f32) {
        ( self.start, self.end, self.value )
    }
}


pub struct BedData {
    pub genome_info: Vec<(String, usize, usize)>, // (chromosome name, length, bin offset)
    pub search: HashMap<String, usize>, // get id for chr
    pub coverage_data: Vec<f32>, // coverage data array for bins
    pub bin_width: f64, // bin width for coverage
    pub threads: usize, // how many worker threads should we use here?
    pub nreads: usize,
}


impl BedData {
    pub fn new( genome_sizes: &[(String, usize)], bin_width:usize ) -> Self {
        let mut offset = 0;
        let mut genome_info=Vec::<(String, usize, usize )>::with_capacity(genome_sizes.len() );
        let mut search = HashMap::new();
        for  (name, size) in  genome_sizes{
            search.insert( name.clone(),  genome_info.len() );
            genome_info.push( (name.clone(), *size, offset + (*size as f64 / bin_width as f64 ).ceil() as usize ) );
            offset += (*size as f64 / bin_width as f64 ).ceil() as usize;
        }
        
        Self{
            genome_info,
            search,
            coverage_data: vec![0.0; genome_sizes.len()],
            bin_width: bin_width as f64,
            threads: 1,
            nreads:0, 
        }
    }
    pub fn add(&mut self, chr:&str, pos:usize, value:f32 ) {
        let id = *self.search.get( chr ).unwrap_or_else( || panic!("chromsome {chr} not found!") );
        self.coverage_data[ self.genome_info[id].2 + (pos as f64 / self.bin_width ) as usize ] += value;
    }

    pub fn write_bigwig( &self, file: &str) -> Result<(),String>{   

        let runtime = tokio::runtime::Builder::new_multi_thread()
             .worker_threads( self.threads )
            .build()
             .expect("Unable to create runtime.");

        let iter = DataIter::new( self );
        let data = BedParserStreamingIterator::wrap_infallible_iter(iter, true);

        let chrom_map: HashMap<String, u32> = self.genome_info.iter().map(|(chrom, len, _)| (chrom.clone(), *len as u32)).collect();

        let outfile = Path::new(file);

        // Create the BigWig writer
        let outb = BigWigWrite::create_file(outfile, chrom_map)
            .map_err(|e| format!("Failed to create BigWig file: {}", e))?;

        // Write data using the runtime
        outb.write(data, runtime)
            .map_err(|e| format!("Failed to write BigWig file: {}", e))?;

        Ok(())
    }

    /// This checks the offsets of all chromosomes and returns None if the id exeeds the range.
    pub fn current_chr_for_id(&self, id:usize) -> Option<( String, usize, usize)> {
        for ( chr, length, offset ) in &self.genome_info{
            if id >= *offset &&  (id - offset) * (self.bin_width as usize) < *length {
                return Some( (chr.clone(), *length, *offset) );
            }
        }
        // the id did not map to any entry here!
        return None
    }
}



pub struct DataIter<'a> {
    /// all data
    pub data: &'a BedData,
    /// the pointer to the BedData::coverage_data id
    pub current_bin: usize,
    /// the current chromosome name and offset for that chromosome
    pub current_chr: Option<( String, usize, usize )>, 
}

impl<'a> DataIter<'a> {
  pub fn new ( data: &'a BedData) ->Self {
    let ret = Some( (
        data.genome_info[0].0.clone(), 
        data.genome_info[0].1,
        data.genome_info[0].2 
      )
    );
    Self{
      data,
      current_bin: 0,
      current_chr: ret,
    }
  }
}



impl<'a> Iterator for DataIter<'a> {
  type Item = (String, bigtools::Value);

  fn next(&mut self) -> Option<Self::Item> {

    let bin_width = self.data.bin_width as usize;
    match &self.current_chr {
      Some((chr, size, offset)) => {

        let mut rel_bin = self.current_bin - offset;
        let start:u32 = (rel_bin * bin_width as usize).try_into().unwrap();
        let mut skipped = false;
        let val = self.data.coverage_data[self.current_bin];
        while self.data.coverage_data[self.current_bin] == val{
          self.current_bin +=1;
          rel_bin +=1;
          skipped = true;
          if (rel_bin * bin_width) >= *size {
            // shit - we overreached out chr!
            // chould not be an issue as we flip back one in the next step
            break;
          }
        }
        if skipped {
          // we need to report the last zero!
          self.current_bin -=1;
        }
       
        rel_bin = self.current_bin - offset;
        // Create the return value based on current values
        let ret = (
          chr.to_string(),
          bigtools::Value {
            start: start,
            end: (rel_bin * bin_width  + bin_width )
              .min(*size)
              .try_into()
              .unwrap(),
            value: val as f32,
          },
        );

        // Increment the bin index
        self.current_bin += 1;

        // Check if the current bin exceeds the size of the chromosome
        let rel_bin = self.current_bin - offset;
        if (rel_bin * bin_width) >= *size {
          // Move to the next chromosome
          self.current_chr = self.data.current_chr_for_id(self.current_bin);
        }
        Some(ret)
        
      }
      None => None, // No more chromosomes to process
    }
  }
}