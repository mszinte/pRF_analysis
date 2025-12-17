#!/bin/bash

# Script to create participants.tsv and participants.json files
# Run from: /scratch/mszinte/data/amblyo7T_prf/
# sh /home/mszinte/projects/pRF_analysis/amblyo7T_prf/preproc/bids/create_participants.sh 

cd /scratch/mszinte/data/amblyo7T_prf/

# Create participants.json
cat > participants.json << 'EOF'
{
	"participant_id": {
		"LongName": "patient identification",
		"Description": "anonymized, BIDS compatible identification"
	},
	"sex": {
		"LongName": "Sex",
		"Description": "Sex of the participant",
		"Levels": { 
			"M": "Male",
			"F": "Female"
		}
	},
	"age": {
		"LongName": "Age",
		"Description": "Age of the subject at time of scan with 3 digits"
	},
	"group": {
		"LongName": "Group",
		"Description": "Participant group classification",
		"Levels": {
			"patient": "Patient with amblyopia",
			"control": "Control participant"
		}
	},
	"type": {
		"LongName": "Amblyopia type",
		"Description": "Type of amblyopia",
		"Levels": {
			"micro": "Microstrabismus",
			"strab": "Strabismus",
			"NA": "Not applicable (control)"
		}
	},
	"amblyopic_eye": {
		"LongName": "Amblyopic eye",
		"Description": "Which eye is amblyopic and type of deviation"
	},
	"alignment_pd": {
		"LongName": "Alignment in prism diopters",
		"Description": "Eye alignment measurement in prism diopters"
	},
	"log_mar": {
		"LongName": "LogMAR visual acuity",
		"Description": "Logarithm of the Minimum Angle of Resolution visual acuity score"
	},
	"patching": {
		"LongName": "Patching history",
		"Description": "History of patching treatment"
	},
	"ops": {
		"LongName": "Operations history",
		"Description": "History of surgical operations"
	},
	"refractive_error_amblyopic": {
		"LongName": "Refractive error (amblyopic eye)",
		"Description": "Refractive error of the amblyopic eye"
	},
	"refractive_error_fellow": {
		"LongName": "Refractive error (fellow eye)",
		"Description": "Refractive error of the fellow eye"
	},
	"log_mar_fellow": {
		"LongName": "LogMAR visual acuity (fellow eye)",
		"Description": "LogMAR visual acuity of the fellow eye"
	}
}
EOF

# Create participants.tsv
cat > participants.tsv << 'EOF'
participant_id	sex	age	group	type	amblyopic_eye	alignment_pd	log_mar	patching	ops	refractive_error_amblyopic	refractive_error_fellow	log_mar_fellow
sub-01	F	38	patient	micro	R Mixed	<4	0.32	No	No	-1.00	-4.00/-0.50x120	0.06
sub-02	F	20	patient	micro	R Micro	<4	0.62	Yes	No	+1.75DS	+1.75DS	0.06
sub-03	M	20	patient	strab	L XOT	6	0.58	Yes	No	+3.50/-175x20	plano	0.14
sub-04	M	67	patient	strab	L SOT	6	0.38	Yes	No	+7.25/-1.75x85	+6.75/-1.75x85	0.08
sub-05	M	38	patient	strab	R SOT	6	0.7	Yes	No	plano	-0.50/-0.50x160	0.00
sub-06	F	19	patient	strab	L SOT	10	0.36	Yes	Yes	+6.25/-0.75x80	+5.25/-0.75x110	0.06
sub-07	M	19	patient	micro	R Micro	<4	0.4	No	No	+4.25/-3.25x17.5	+0.50/-0.50x170	0.08
sub-08	M	24	patient	strab	L XOT	8-14	1	Yes	No	+4.00/-1.50x180	-1.50DS	0.10
sub-09	M	25	patient	strab	R SOT	9	0.72	Yes	No	+1.25/-0.75x165	-3.25/-1.00x180	0.02
sub-10	F	48	patient	strab	R SOT	12	0.56	Yes	yes	plano/-0.75x80	NA	NA
sub-11	F	22	patient	strab	L XOT	6	0.42	Yes	yes	+0.25DS	+0.50/-0.25x110	0.06
sub-12	M	48	patient	micro	R mixed	<4	0.38	Yes	No	+7.00/-2.00x10	-3.50/-0.75x145	0.08
EOF

echo "participants.json and participants.tsv created successfully!"
