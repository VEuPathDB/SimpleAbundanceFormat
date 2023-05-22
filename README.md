# SimpleAbundanceFormat
Formerly known as PopBio-interchange-format. This is a new repo that will build on https://github.com/chowington/PopBio-interchange-format and add 2.0 features.

This is VEuPathDB's implementation of MIReAD (https://www.nature.com/articles/s41597-019-0042-5) as the prefered file format for abundance and pathogen surveillance data. Internally we convert it to ISA-Tab with the "SAF-Wizard" and then load into our system.

Description of Simple Abundance Format (SAF)
--------------------------------------------
Field Name |Format|Requirement|Details|Multi-valued?
-----------|------|-------|-----------|-------
collection_ID|string|Mandatory|Identifier for collection event e.g. ABC_2018_collection_00001
sample_ID|string|Mandatory|Identifier for sample event e.g. ABC_2018_sample_00001
collection_start_date|ISO 8601 date format (YYYY-MM-DD)|Mandatory|Date at which traps deployed
collection_end_date|ISO 8601 date format (YYYY-MM-DD)|Mandatory|Date at which traps collected and specimens removed for processing
trap_ID|string|Advisory|Internal trap ID for the collection, may be used as part of processing for distinct collection events
GPS_latitude|GPS decimal degrees|Mandatory|Latitude for collection site, max. 6 decimal places
GPS_longitude|GPS decimal degrees|Mandatory|Longitude for collection site, max. 6 decimal places
GPS_qualifier|string|Optional|controlled vocabulary/ontology term to describe source and precision of GPS coordinates
location_description|string|Optional|Collection location description e.g. Orlando
location_ADM2|string|Optional|Administrative level 2 for collection e.g. Orange County
location_ADM1|string|Optional|Administrative level 1 for collection e.g. Florida
location_country|string|Optional|Country in which collection occurred e.g. United States of America
trap_type|string|Mandatory|Trap type e.g. CDC light trap, New Jersey Trap, defined further in protocols section
attractant|string|Advisory|List of attractants used in the trap e.g. CO2, light|yes, semicolon delimited
trap_number|integer|Mandatory|Number of traps deployed (Default is 1)
trap_duration|integer|Mandatory|Number of nights/days trap was deployed (Default is 1)
species|string|Species|binomial species name for collected specimens
species_identification_method|string|Mandatory|Protocol for asserting species identification
developmental_stage|string|Mandatory|developmental stage e.g. adult
combined_feeding_and_gonotrophic_status|string|Optional|e.g. fed, unfed, semi-gravid, gravid
sex|string|Mandatory|sex of specimens e.g. female/male/mixed
sample_count|integer|Mandatory|count of specimens from collection
collection_comment|string|Optional|free text comment about the collection event
sample_comment|string|Optional|free text comment about the sample material
species_comment|string|Optional|free text comments about the species identification process
phenotypes|string|Optional| phenotype should be formatted as follows and separated by vertical bars '\|' <br/><br/>**PCR_VIVAX;arthropod infection status;Plasmodium vivax;present**<br/><br/>Spaces before/after the semicolon and vertical bar delimiters are allowed.<br/><br/>PCR_VIVAX is the protocol type (defined in study_protocols in config file). The next three values are the Observable, Attribute and Value that describe the phenotype (GMOD Chado style) and each value must be a term described in the study_terms section of the config file. Currently only pathogen infection status phenotypes are supported. The value for positive assays must match /^(?:present\|positive\|confirmed\|detected)/i|yes, pipe delimited

SAF2.0 additions

Field Name |Format|Requirement|Details
-----------|------|-------|-----------
location_comment|string|Optional|free text comment about the collection site
any_other_column_name|string or number|Optional|must be specified in config file, see below



Example configuration file
--------------------------
Explanations in comments
```
study_identifier : 2023-maine-ricinus
study_title : Tick surveillance from Maine
study_submission_date : 2023
study_description : Tick surveillance from the University of Maine Cooperative Extension, Tick Lab. Residents of Maine may submit ticks for species IDs of the ticks, and pathogen testing. 

# This section, and the study_contacts sections are actually
# just YAML representations of the ISA-Tab investigation sheet sections.
# Additional publications/contacts can be added with indented '-' sections (see study_protocols)
study_publications :
  - study_publication_doi : 
    study_pubmed_id : 
    study_publication_author_list :  University of Maine Cooperative Extension, Tick Lab
    study_publication_title : Tick submission reports
    study_publication_status : unpublished
    study_publication_status_term_source_ref :
    study_publication_status_term_accession_number :

study_contacts :
  - study_person_last_name : Smith
    study_person_first_name : John
    study_person_email : John Smith <not.real.email@maine.edu>
    study_person_affiliation : University of Maine, Cooperative Extension Tick Lab

# these are the species that will be found in the main species column
study_species :
  Ixodidae : NCBITaxon_6939
  Ixodes scapularis : NCBITaxon_6945
  Dermacentor variabilis : NCBITaxon_34621  

# and the developmental_stage column
study_developmental_stages :
  Nymph : UBERON_0014405
  Adult : UBERON_0018241
  Larvae : UBERON_0000069
  N/A : OBI_0000852
  Unknown : OBI_0000852

# and the sex column
study_sexes :
  Male : PATO_0000384
  Female : PATO_0000383
  N/A : OBI_0000852
  Unknown : OBI_0000852

# These protocols are also already in ISA-Tab format.
# You put the study_protocol_name in the species_identification_method, trap_type and phenotype columns
study_protocols :
  - study_protocol_name : CITIZEN
    study_protocol_type : citizen science collection
    study_protocol_type_term_source_ref : EUPATH
    study_protocol_type_term_accession_number : EUPATH_0043237
    study_protocol_description :  Residents of the state of Maine, USA, sent in ticks they find (usually on themself) to a central laboratory for analysis.

  - study_protocol_name : SPECIES
    study_protocol_type : species identification method
    study_protocol_type_term_source_ref : EUPATH
    study_protocol_type_term_accession_number : OBI_0001624
    study_protocol_description : Detailed species identification methods were not described.

  - study_protocol_name : PCR
    study_protocol_type : real time PCR
    study_protocol_type_term_source_ref : EUPATH
    study_protocol_type_term_accession_number : OBI_0001624
    study_protocol_description : Screening of tick-borne pathogens was performed using PCR.

# any additional values used in the input data need to be mapped to terms here
study_terms :
  pool : OBI_0302716
  Present : PATO_0000467
  Absent : PATO_0000462
  Borrelia burgdorferi sensu lato : NCBITaxon_139
  Borrelia miyamotoi : NCBITaxon_47466
  Babesia microti : NCBITaxon_5868

# SAF2.0 addition
columns :
  # mark a required column as not-required (not recommended!)
  # with an optional default value to be used instead
  - heading : trap_duration
    required : false
    default : 1

  # a totally new column (will be required: true by default)
  # values must have ontology terms provided in study_terms when allow_freetext is false
  - heading : host_species
    required : true
    type : string
    describes : collection
    controlled : false
```

Custom column configuration fields

Field | Default | Details
------|---------|--------
heading | none | the column heading, exactly as it appears; you may add new columns this way
required | true | true/false; use this to disable mandatory columns (not recommended!)
default | none | if the column is missing or if there is no value in a particular row, this value will be used instead
type | string | string/number/date; type of value expected (this will be validated as far as possible)
describes | sample | which aspect of the data does this column describe? must be exactly one of these: <br/> `location`, `collection`, `sample`, `organism identification assay`, `genotyping assay`, `insecticide resistance assay`, `pathogen detection assay`, `blood meal assay`
controlled | true | true/false; for string type columns, must there be ontology terms in study_terms for all values?
protocol | none | only applies to (and is required for) `describes: xxx assay` columns, use the study_protocol_name from the study_protocols section (adding a section if required), e.g. `PCR` in the above example  

See default-column-config.yaml for the built-in column definitions
